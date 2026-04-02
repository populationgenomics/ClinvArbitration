"""
Publish a new versioned record on Zenodo. Mostly generated via Claude.

Creates a new version of an existing Zenodo record, replaces its files with a
local tarball, updates the version label, removes  and publishes.

Usage:
    python publish_to_zenodo.py <record_id> <gcp_secret_name> <file_path>

Arguments:
    record    Integer ID of an existing version of the Zenodo record
    secret    Name of the GCP Secret Manager secret holding the Zenodo token
    tarball   Path to the file to attach to the new record version
    success   Path to an output file, to contain the new record URL
"""

import argparse
import sys
import zoneinfo
from datetime import datetime
from pathlib import Path

import google.auth
import requests
from google.cloud import secretmanager

TIMEZONE = zoneinfo.ZoneInfo('Australia/Brisbane')
TODAY = datetime.now(tz=TIMEZONE)
TODAY_YM = TODAY.strftime('%B %Y')
TODAY_YMD = TODAY.strftime('%Y-%m-%d')
ZENODO_BASE = 'https://zenodo.org/api'


def get_secret(secret_name: str) -> str:
    """Retrieve a secret value from GCP Secret Manager using ADC."""
    _, project_id = google.auth.default()
    if not project_id:
        raise RuntimeError(
            'Could not determine GCP project ID from application default credentials. '
            'Run `gcloud auth application-default login` or set GOOGLE_CLOUD_PROJECT.',
        )
    client = secretmanager.SecretManagerServiceClient()
    secret_path = f'projects/{project_id}/secrets/{secret_name}/versions/latest'
    response = client.access_secret_version(name=secret_path)
    return response.payload.data.decode('utf-8').strip()


def _check(response: requests.Response, action: str) -> dict:
    """Raise with context on non-2xx responses."""
    if not response.ok:
        raise RuntimeError(
            f'{action} failed [{response.status_code}]: {response.text[:500]}',
        )
    return response.json()


def get_deposition(record_id: int, token: str) -> dict:
    """Fetch an existing deposition to confirm it exists."""
    print(f'Fetching existing record {record_id}...')
    r = requests.get(
        f'{ZENODO_BASE}/deposit/depositions/{record_id}',
        params={'access_token': token},
        timeout=30,
    )
    return _check(r, f'GET deposition {record_id}')


def create_new_version(record_id: int, token: str) -> dict:
    """Create a new draft version of an existing deposition."""
    print('Creating new version draft...')
    r = requests.post(
        f'{ZENODO_BASE}/deposit/depositions/{record_id}/actions/newversion',
        params={'access_token': token},
        timeout=30,
    )
    return _check(r, 'Create new version')


def get_draft_id(new_version_response: dict) -> int:
    """Extract the new draft deposition ID from the latest_draft URL."""
    latest_draft_url: str = new_version_response['links']['latest_draft']
    return int(latest_draft_url.rstrip('/').split('/')[-1])


def clear_inherited_files(draft_id: int, token: str) -> None:
    """Delete all files inherited from the previous version."""
    print('Clearing inherited files...')
    r = requests.get(
        f'{ZENODO_BASE}/deposit/depositions/{draft_id}/files',
        params={'access_token': token},
        timeout=30,
    )
    files = _check(r, 'List draft files')
    for f in files:
        file_id = f['id']
        rd = requests.delete(
            f'{ZENODO_BASE}/deposit/depositions/{draft_id}/files/{file_id}',
            params={'access_token': token},
            timeout=30,
        )
        if not rd.ok:
            raise RuntimeError(
                f'Delete file {file_id} failed [{rd.status_code}]: {rd.text[:500]}',
            )
        print(f'  Deleted inherited file: {f["filename"]}')


def upload_file(draft: dict, file_path: Path, token: str) -> None:
    """Upload a file to the draft using the bucket API."""
    bucket_url: str = draft['links']['bucket']
    filename = file_path.name
    print(f'Uploading {filename} ({file_path.stat().st_size / 1024 / 1024:.1f} MB)...')
    with file_path.open('rb') as fh:
        r = requests.put(
            f'{bucket_url}/{filename}',
            params={'access_token': token},
            data=fh,
            timeout=600,
        )
    _check(r, f'Upload {filename}')
    print('  Upload complete.')


def update_metadata(draft_id: int, existing_metadata: dict, token: str) -> None:
    """Update only the version field in the record metadata."""
    version_label = f"ClinvArbitration data release - {TODAY_YM}"
    print(f'Setting version label: "{version_label}"')
    metadata = dict(existing_metadata)
    metadata |= {
        'version': version_label,
        'date': TODAY_YMD,
        'publication_date': TODAY_YMD,
    }
    _ = metadata.pop('dates', None)
    _ = metadata.pop('doi', None)
    _ = metadata.pop('prereserve_doi', None)
    r = requests.put(
        f'{ZENODO_BASE}/deposit/depositions/{draft_id}',
        params={'access_token': token},
        json={'metadata': metadata},
        timeout=30,
    )
    _check(r, 'Update metadata')


def publish(draft_id: int, token: str) -> dict:
    """Publish the draft, making it a live record."""
    print('Publishing...')
    r = requests.post(
        f'{ZENODO_BASE}/deposit/depositions/{draft_id}/actions/publish',
        params={'access_token': token},
        timeout=60,
    )
    return _check(r, 'Publish')


def main(record: int, token_secret: str, tarball_path: str, success_file: str) -> None:
    tarball = Path(tarball_path)
    if not tarball.exists():
        raise FileNotFoundError(f'File not found: {tarball}')

    print(f'Retrieving Zenodo token from GCP secret "{token_secret}"...')
    token = get_secret(token_secret)

    existing = get_deposition(record, token)
    new_version_resp = create_new_version(record, token)
    draft_id = get_draft_id(new_version_resp)
    print(f'New draft deposition ID: {draft_id}')

    # Fetch the full draft to get the bucket URL and current metadata
    draft = get_deposition(draft_id, token)

    clear_inherited_files(draft_id, token)
    upload_file(draft, tarball, token)
    update_metadata(draft_id, existing['metadata'], token)
    result = publish(draft_id, token)

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
    parser.add_argument('--secret', help='GCP name for the Zenodo token secret', required=True)
    parser.add_argument('--tarball', type=str, help='path to the file to upload', required=True)
    parser.add_argument('--success', type=str, help='if successful, write new record URL', required=True)
    args = parser.parse_args()

    try:
        main(args.record, args.secret, args.tarball, args.success)
    except Exception as exc:  # noqa: BLE001
        print(f'Error: {exc}', file=sys.stderr)
        sys.exit(1)
