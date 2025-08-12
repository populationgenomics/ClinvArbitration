from typing import TYPE_CHECKING

from cpg_utils import Path, config, hail_batch

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def package_data_for_release(
    pm5: dict[str, Path],
    clinvar_decisions: dict[str, Path],
    output: Path,
) -> 'BashJob':
    """Localise all the previously generated data into a folder - tarball it, and write out as a single file."""

    batch_instance = hail_batch.get_batch('Run ClinvArbitration')

    job = batch_instance.new_bash_job('PackageForRelease')
    job.image(config.config_retrieve(['workflow', 'driver_image'])).storage('10G')

    decisions_ht = batch_instance.read_input(clinvar_decisions['clinvar_decisions'])
    decisions_tsv = batch_instance.read_input(clinvar_decisions['tsv'])
    pm5_ht = batch_instance.read_input(pm5['ht'])
    pm5_tsv = batch_instance.read_input(pm5['tsv'])

    # create a new folder, move the files into it
    # unpack all the already compressed HTs into it
    # then compress once, containing all files for distribution
    job.command(
        f"""
        mkdir clinvarbitration_data
        tar -xf {pm5_ht} -C clinvarbitration_data
        tar -xf {decisions_ht} -C clinvarbitration_data
        mv {pm5_tsv} clinvarbitration_data/clinvar_decisions.pm5.tsv
        mv {decisions_tsv} clinvarbitration_data/clinvar_decisions.tsv
        tar -czf {job.output} \
            clinvarbitration_data/clinvar_decisions.ht \
            clinvarbitration_data/clinvar_decisions.pm5.ht \
            clinvarbitration_data/clinvar_decisions.pm5.tsv \
            clinvarbitration_data/clinvar_decisions.tsv
    """,
    )
    batch_instance.write_output(job.output, output)

    return job
