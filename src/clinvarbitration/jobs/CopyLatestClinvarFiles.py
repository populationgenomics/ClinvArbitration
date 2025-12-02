from typing import TYPE_CHECKING

from cpg_utils import Path

from clinvarbitration.cpg_internal.utils import make_me_a_job

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

    job = make_me_a_job('CopyLatestClinvarFiles').storage('10Gi')

    directory = 'https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/'

    for filename, output_name in [('submission_summary.txt.gz', submissions), ('variant_summary.txt.gz', variants)]:

        job.command(f'wget -q {directory}{filename} -O - | gcloud storage cp - {output_name}')

    return job
