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
    bash_job = get_batch().new_bash_job('CopyLatestClinvarFiles')
    bash_job.image(config_retrieve(['workflow', 'driver_image']))

    directory = 'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/'
    sub_file = 'submission_summary.txt.gz'
    var_file = 'variant_summary.txt.gz'

    bash_job.command(f'wget -q {directory}{sub_file} -O {bash_job.subs}')
    bash_job.command(f'wget -q {directory}{var_file} -O {bash_job.vars}')

    get_batch().write_output(bash_job.subs, submissions)
    get_batch().write_output(bash_job.vars, variants)

    return bash_job
