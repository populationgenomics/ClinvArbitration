from typing import TYPE_CHECKING

from cpg_utils import Path, config, hail_batch

from clinvarbitration.cpg_internal.utils import make_me_a_job

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def create_new_release(
    tarball: Path,
    successfile: Path,
) -> 'BashJob':
    """Localise the Tarball, create and publish a new release, write success file."""

    batch_instance = hail_batch.get_batch('Run ClinvArbitration')

    job = make_me_a_job('PublishToZenodo').storage('10G')

    tarball_local = batch_instance.read_input(tarball)

    # a record ID to create a new version of, not necessarily the latest
    zenodo_record = config.config_retrieve(['workflow', 'zenodo_id'])

    # name of the GCP ID to extract a zenodo auth token from
    zenodo_secret = config.config_retrieve(['workflow', 'zenodo_secret'])

    job.command(
        f"""
        set -euo pipefail
        python -m clinvarbitration.scripts.publish_to_zenodo \\
            --record {zenodo_record} \\
            --secret {zenodo_secret} \\
            --tarball {tarball_local} \\
            --success {job.output}
    """,
    )

    batch_instance.write_output(job.output, successfile)

    return job
