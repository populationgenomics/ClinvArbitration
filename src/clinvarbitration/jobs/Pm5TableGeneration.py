from typing import TYPE_CHECKING

from cpg_utils import Path
from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import get_batch

from clinvarbitration.scripts import clinvar_by_codon

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def generate_pm5_data(
    annotated_snvs: Path,
    output_root: Path,
) -> 'BashJob':
    """Generate PM5 data (index pathogenic missense variants by codon/transcript)."""

    annotated_snvs_local = get_batch().read_input(annotated_snvs)

    job = get_batch().new_bash_job('Pm5TableGeneration')
    job.image(config_retrieve(['workflow', 'driver_image'])).storage('10G')

    job.declare_resource_group(output={'ht.tar': '{root}.ht.tar', 'json': '{root}.json'})

    # write both HT and JSON outputs to the same root location
    job.command(f'python3 {clinvar_by_codon.__file__} -i {annotated_snvs_local} -o {job.output}')

    # compress the HT and remove as a single file
    job.command(
        f'mv {job.output}.ht clinvar_decisions.pm5.ht && tar -cf {job.output}.ht.tar clinvar_decisions.pm5.ht',
    )

    # write both outputs together
    get_batch().write_output(job.output, output_root)

    return job
