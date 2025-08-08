from typing import TYPE_CHECKING

from cpg_utils import Path, config, hail_batch

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def copy_latest_files(
    submissions: Path,
    variants: Path,
) -> 'BashJob':
    """
    gets the remote resources for submissions and variants

    Args:
        submissions (Pathlike): path to write the submission file
        variants (Pathlike): path to write the variant file
    """
    batch_instance = hail_batch.get_batch('Run ClinvArbitration')

    job = batch_instance.new_bash_job('CopyLatestClinvarFiles')

    job.image(config.config_retrieve(['workflow', 'driver_image']))

    directory = 'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/'
    sub_file = 'submission_summary.txt.gz'
    var_file = 'variant_summary.txt.gz'

    job.command(f'wget -q {directory}{sub_file} -O {job.subs}')
    job.command(f'wget -q {directory}{var_file} -O {job.vars}')

    batch_instance.write_output(job.subs, submissions)
    batch_instance.write_output(job.vars, variants)

    return job
