"""
Publish a new versioned record on Zenodo. Mostly generated via Claude.

Creates a new version of an existing Zenodo record (or resumes an unpublished draft
left by a previous failed run), replaces its files with a local tarball, updates the
version label, and publishes.

Usage:
    python publish_to_zenodo.py <record_id> <gcp_secret_name> <file_path>

Arguments:
    record    Integer ID of an existing version of the Zenodo record
    secret    Name of the GCP Secret Manager secret holding the Zenodo token
    project   ID of the GCP project ID holding the Zenodo secret
    tarball   Path to the file to attach to the new record version
    success   Path to an output file, to contain the new record URL
    token     A naked auth token, alternative to secret + project. Use caution in a logged environment
"""

import argparse
import hashlib
import sys
import time
import zoneinfo
from datetime import datetime
from pathlib import Path

import requests
from google.cloud import secretmanager
from requests.adapters import HTTPAdapter, Retry

TIMEZONE = zoneinfo.ZoneInfo('Australia/Brisbane')
TODAY = datetime.now(tz=TIMEZONE)
TODAY_YM = TODAY.strftime('%B %Y')
TODAY_YMD = TODAY.strftime('%Y-%m-%d')
ZENODO_BASE = 'https://zenodo.org/api'

# Zenodo intermittently drops long-running uploads (seen as SSLEOFError / ConnectionResetError),
# so the file upload is retried explicitly, and every attempt is verified against a local MD5.
UPLOAD_ATTEMPTS = 5
UPLOAD_TIMEOUT = 600
RETRYABLE_STATUSES = (429, 500, 502, 503, 504)
RETRYABLE_ERRORS = (requests.ConnectionError, requests.Timeout, requests.exceptions.ChunkedEncodingError)
MD5_CHUNK = 8 * 1024 * 1024


def make_session() -> requests.Session:
    """
    A session with automatic retries for the small metadata calls.

    Only idempotent methods are retried after a read error; POST (newversion, publish) is retried
    on connection failures alone, so a lost response can't create a duplicate draft or publish twice.
    The large file upload deliberately bypasses this and uses its own loop in upload_file().
    """
    retry = Retry(
        total=5,
        connect=5,
        read=3,
        status=5,
        backoff_factor=2,
        status_forcelist=RETRYABLE_STATUSES,
        allowed_methods=frozenset({'GET', 'PUT', 'DELETE'}),
        raise_on_status=False,
    )
    session = requests.Session()
    session.mount('https://', HTTPAdapter(max_retries=retry))
    return session


def md5_of(file_path: Path) -> str:
    """Hex MD5 of a file, read in chunks."""
    digest = hashlib.md5(usedforsecurity=False)
    with file_path.open('rb') as fh:
        for chunk in iter(lambda: fh.read(MD5_CHUNK), b''):
            digest.update(chunk)
    return digest.hexdigest()


def get_secret(secret_name: str, secret_project: str) -> str:
    """Retrieve a secret value from GCP Secret Manager using ADC."""
    client = secretmanager.SecretManagerServiceClient()
    secret_path = f'projects/{secret_project}/secrets/{secret_name}/versions/latest'
    response = client.access_secret_version(name=secret_path)
    return response.payload.data.decode('utf-8').strip()


def _check(response: requests.Response, action: str) -> dict:
    """Raise with context on non-2xx responses."""
    if not response.ok:
        raise RuntimeError(
            f'{action} failed [{response.status_code}]: {response.text[:500]}',
        )
    return response.json()


def get_deposition(session: requests.Session, record_id: int, token: str) -> dict:
    """Fetch an existing deposition to confirm it exists."""
    print(f'Fetching existing record {record_id}...')
    r = session.get(
        f'{ZENODO_BASE}/deposit/depositions/{record_id}',
        params={'access_token': token},
        timeout=30,
    )
    return _check(r, f'GET deposition {record_id}')


def create_new_version(session: requests.Session, record_id: int, token: str) -> dict:
    """Create a new draft version of an existing deposition."""
    print('Creating new version draft...')
    r = session.post(
        f'{ZENODO_BASE}/deposit/depositions/{record_id}/actions/newversion',
        params={'access_token': token},
        timeout=30,
    )
    return _check(r, 'Create new version')


def get_draft_id(deposition: dict) -> int | None:
    """
    Extract the draft deposition ID from a deposition's latest_draft link.

    Zenodo includes this link on both the newversion response and on a published record which
    already has an unpublished draft. Returns None when no such link is present.
    """
    latest_draft_url: str | None = deposition.get('links', {}).get('latest_draft')
    if not latest_draft_url:
        return None
    return int(latest_draft_url.rstrip('/').split('/')[-1])


def find_existing_draft(session: requests.Session, existing: dict, token: str) -> int | None:
    """
    Return the ID of an unpublished draft already attached to this record, if there is one.

    Lets a failed run (e.g. a dropped upload) be resumed against the same draft rather than
    leaving an orphaned draft behind and attempting to create a second one.
    """
    draft_id = get_draft_id(existing)
    if draft_id is None:
        return None
    draft = get_deposition(session, draft_id, token)
    if draft.get('submitted', False):
        # link points at something already published, so it isn't a resumable draft
        return None
    return draft_id


def get_or_create_draft(session: requests.Session, record_id: int, existing: dict, token: str) -> int:
    """Reuse an existing unpublished draft of this record, or create a new version draft."""
    draft_id = find_existing_draft(session, existing, token)
    if draft_id is not None:
        print(f'Found existing unpublished draft {draft_id}, resuming with it.')
        return draft_id

    new_version_resp = create_new_version(session, record_id, token)
    draft_id = get_draft_id(new_version_resp)
    if draft_id is None:
        raise RuntimeError(f'New version response had no latest_draft link: {new_version_resp.get("links")}')
    print(f'New draft deposition ID: {draft_id}')
    return draft_id


def clear_inherited_files(session: requests.Session, draft_id: int, token: str) -> None:
    """Delete all files inherited from the previous version."""
    print('Clearing inherited files...')
    r = session.get(
        f'{ZENODO_BASE}/deposit/depositions/{draft_id}/files',
        params={'access_token': token},
        timeout=30,
    )
    files = _check(r, 'List draft files')
    for f in files:
        file_id = f['id']
        rd = session.delete(
            f'{ZENODO_BASE}/deposit/depositions/{draft_id}/files/{file_id}',
            params={'access_token': token},
            timeout=30,
        )
        if not rd.ok:
            raise RuntimeError(
                f'Delete file {file_id} failed [{rd.status_code}]: {rd.text[:500]}',
            )
        print(f'  Deleted inherited file: {f["filename"]}')


def _attempt_upload(url: str, file_path: Path, token: str, local_md5: str) -> str | None:
    """
    One upload attempt. Returns None on verified success, or a description of the problem to retry.

    Raises RuntimeError on non-retryable HTTP errors (4xx), since those won't improve on retry.
    Uses a bare requests.put rather than the retrying session so attempts aren't multiplied.
    """
    try:
        with file_path.open('rb') as fh:
            r = requests.put(url, params={'access_token': token}, data=fh, timeout=UPLOAD_TIMEOUT)
    except RETRYABLE_ERRORS as exc:
        return f'connection error: {exc}'

    if r.status_code in RETRYABLE_STATUSES:
        return f'HTTP {r.status_code}: {r.text[:200]}'

    remote_md5 = _check(r, f'Upload {file_path.name}').get('checksum', '').removeprefix('md5:')
    if remote_md5 != local_md5:
        return f'checksum mismatch (local {local_md5}, Zenodo {remote_md5 or "absent"})'
    return None


def upload_file(draft: dict, file_path: Path, token: str) -> None:
    """
    Upload a file to the draft using the bucket API, retrying on dropped connections.

    The bucket API overwrites by filename, so a retry after a partial upload is safe. Each
    attempt reopens the file from the start and is only accepted if Zenodo's reported MD5
    matches the local file.
    """
    url = f'{draft["links"]["bucket"]}/{file_path.name}'
    size_mb = file_path.stat().st_size / 1024 / 1024

    print(f'Computing MD5 of {file_path.name}...')
    local_md5 = md5_of(file_path)

    for attempt in range(1, UPLOAD_ATTEMPTS + 1):
        print(f'Uploading {file_path.name} ({size_mb:.1f} MB), attempt {attempt}/{UPLOAD_ATTEMPTS}...')
        problem = _attempt_upload(url, file_path, token, local_md5)
        if problem is None:
            print(f'  Upload complete, checksum verified ({local_md5}).')
            return
        if attempt == UPLOAD_ATTEMPTS:
            raise RuntimeError(f'Upload {file_path.name} failed after {UPLOAD_ATTEMPTS} attempts; last: {problem}')
        delay = 5 * 2**attempt
        print(f'  Attempt {attempt} failed ({problem}); retrying in {delay}s...')
        time.sleep(delay)


def update_metadata(session: requests.Session, draft_id: int, existing_metadata: dict, token: str) -> None:
    """Update only the version field in the record metadata."""
    version_label = f'ClinvArbitration data release - {TODAY_YM}'
    print(f'Setting version label: "{version_label}"')
    metadata = dict(existing_metadata)
    metadata |= {
        'version': version_label,
        'title': version_label,
        'date': TODAY_YMD,
        'publication_date': TODAY_YMD,
    }
    _ = metadata.pop('dates', None)
    _ = metadata.pop('doi', None)
    _ = metadata.pop('prereserve_doi', None)
    r = session.put(
        f'{ZENODO_BASE}/deposit/depositions/{draft_id}',
        params={'access_token': token},
        json={'metadata': metadata},
        timeout=30,
    )
    _check(r, 'Update metadata')


def publish(session: requests.Session, draft_id: int, token: str) -> dict:
    """Publish the draft, making it a live record."""
    print('Publishing...')
    r = session.post(
        f'{ZENODO_BASE}/deposit/depositions/{draft_id}/actions/publish',
        params={'access_token': token},
        timeout=60,
    )
    return _check(r, 'Publish')


def main(
    record: int,
    token_secret: str,
    token_project: str,
    tarball_path: str,
    success_file: str,
    the_token: str | None = None,
) -> None:
    tarball = Path(tarball_path)
    if not tarball.exists():
        raise FileNotFoundError(f'File not found: {tarball}')

    if the_token:
        print('Using provided token')
        token = the_token
    else:
        print(f'Retrieving Zenodo token from GCP secret "{token_secret}"...')
        token = get_secret(token_secret, secret_project=token_project)

    session = make_session()

    existing = get_deposition(session, record, token)
    draft_id = get_or_create_draft(session, record, existing, token)

    # Fetch the full draft to get the bucket URL and current metadata
    draft = get_deposition(session, draft_id, token)

    # Every step from here is safe to repeat on a resumed draft: listing/deleting files is a no-op
    # once cleared, the bucket upload overwrites by filename, and the metadata PUT is idempotent.
    clear_inherited_files(session, draft_id, token)
    upload_file(draft, tarball, token)
    update_metadata(session, draft_id, existing['metadata'], token)
    result = publish(session, draft_id, token)

    doi = result.get('doi', 'unknown')
    record_url = result.get('links', {}).get('record_html', 'unknown')
    print('\nPublished successfully.')
    print(f'  DOI:  {doi}')
    print(f'  URL:  {record_url}')

    with open(success_file, 'w') as f:
        f.write(record_url)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Publish a new versioned Zenodo record from a local file.',
    )
    parser.add_argument('--record', type=int, help='ID of the existing Zenodo record', required=True)
    parser.add_argument('--secret', help='GCP name for the Zenodo token secret')
    parser.add_argument('--project', help='GCP ID for the Zenodo token project')
    parser.add_argument('--tarball', type=str, help='path to the file to upload', required=True)
    parser.add_argument('--success', type=str, help='if successful, write new record URL', required=True)
    parser.add_argument('--token', type=str, help='the full token, avoiding secret retrieval')
    args = parser.parse_args()

    # if secret and project are provided, dismiss token (probably already None)
    if args.secret and args.project:
        args.token = None

    try:
        main(
            record=args.record,
            token_secret=args.secret,
            token_project=args.project,
            the_token=args.token,
            tarball_path=args.tarball,
            success_file=args.success,
        )
    except Exception as exc:  # noqa: BLE001
        print(f'Error: {exc}', file=sys.stderr)
        sys.exit(1)
