from typing import TYPE_CHECKING

from cpg_utils import config, hail_batch


if TYPE_CHECKING:
    from hailtop.batch.job import BashJob

def make_me_a_job(name: str, attributes: dict[str, str] | None = None) -> 'BashJob':
    """Creates a job using the driver image and inserting a gcloud instruction."""
    job = hail_batch.get_batch().new_bash_job(name=name, attributes=attributes)
    job.image(config.config_retrieve(['workflow', 'driver_image']))
    job.command('gcloud config set storage/parallel_composite_upload_enabled False')
    return job
