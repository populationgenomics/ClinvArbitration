from typing import TYPE_CHECKING

from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import get_batch


if TYPE_CHECKING:
    from hailtop.batch.job import BashJob
    from cpg_utils import Path


def copy_latest_files(
    submissions: 'Path',
    variants: 'Path',
) -> 'BashJob':
    """
    gets the remote resources for submissions and variants

    Args:
        submissions (Pathlike): path to write the submission file
        variants (Pathlike): path to write the variant file
    """
    job = get_batch().new_bash_job('CopyLatestClinvarFiles')
    job.image(config_retrieve(['workflow', 'driver_image']))

    directory = 'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/'
    sub_file = 'submission_summary.txt.gz'
    var_file = 'variant_summary.txt.gz'

    job.command(f'wget -q {directory}{sub_file} -O {job.subs}')
    job.command(f'wget -q {directory}{var_file} -O {job.vars}')

    get_batch().write_output(job.subs, submissions)
    get_batch().write_output(job.vars, variants)

    return job
