from typing import TYPE_CHECKING

from cpg_utils import Path, config, hail_batch

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def generate_pm5_data(
    annotated_snvs: Path,
    output_root: Path,
) -> 'BashJob':
    """Generate PM5 data (index pathogenic missense variants by codon/transcript)."""

    batch_instance = hail_batch.get_batch('Run ClinvArbitration')

    annotated_snvs_local = batch_instance.read_input(annotated_snvs)

    job = batch_instance.new_bash_job('Pm5TableGeneration')
    job.image(config.config_retrieve(['workflow', 'driver_image'])).storage('10G')

    job.declare_resource_group(output={'ht.tar': '{root}.ht.tar', 'json': '{root}.json'})

    # write both HT and JSON outputs to the same root location
    job.command(f'python3 -m clinvarbitration.scripts.clinvar_by_codon -i {annotated_snvs_local} -o {job.output}')

    # compress the HT and remove as a single file
    job.command(
        f'mv {job.output}.ht clinvar_decisions.pm5.ht && tar -cf {job.output}.ht.tar clinvar_decisions.pm5.ht',
    )

    # write both outputs together
    batch_instance.write_output(job.output, output_root)

    return job
