from typing import TYPE_CHECKING

from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import get_batch

from clinvarbitration.scripts import resummarise_clinvar


if TYPE_CHECKING:
    from hailtop.batch.job import BashJob
    from cpg_utils import Path


def generate_new_summary(
    var_file: 'Path',
    sub_file: 'Path',
    output_root: 'Path',
) -> 'BashJob':
    """Gets the remote resources for submissions and variants."""

    job = get_batch().new_bash_job('GenerateNewClinvarSummary')
    job.image(config_retrieve(['workflow', 'driver_image'])).memory('highmem').cpu('2')

    if sites_to_blacklist := config_retrieve(['workflow', 'site_blacklist'], []):
        blacklist_sites = ' '.join(f'"{site}"' for site in sites_to_blacklist)
        blacklist_string = f' -b {blacklist_sites}'
    else:
        blacklist_string = ''

    var_file_local = get_batch().read_input(var_file)
    sub_file_local = get_batch().read_input(sub_file)

    job.declare_resource_group(
        output={
            'ht.tar': '{root}.ht.tar',
            'vcf.bgz': '{root}.vcf.bgz',
            'vcf.bgz.tbi': '{root}.vcf.bgz.tbi',
        },
    )

    job.command(
        f'python3 {resummarise_clinvar.__file__} \
        -v {var_file_local} \
        -s {sub_file_local} \
        -o {job.output} {blacklist_string}',
    )

    # don't tar from current location, we'll catch all the tmp pathing
    job.command(f'mv {job.output}.ht clinvar_decisions.ht && tar -cf {job.output}.ht.tar clinvar_decisions.ht')

    # selectively copy back some outputs
    get_batch().write_output(job.output, output_root)

    return job
