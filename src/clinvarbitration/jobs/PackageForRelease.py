from typing import TYPE_CHECKING

from cpg_utils import Path, config, hail_batch

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def package_data_for_release(
    pm5_json: Path,
    pm5_ht: Path,
    clinvar_decisions: Path,
    output: Path,
) -> 'BashJob':
    """
    Localise all the previously generated data into a folder
    tarball it, and write out as a single file
    """

    batch_instance = hail_batch.get_batch('Run ClinvArbitration')

    job = batch_instance.new_bash_job('PackageForRelease')
    job.image(config.config_retrieve(['workflow', 'driver_image'])).storage('10G')

    pm5_ht_local = batch_instance.read_input(pm5_ht)
    pm5_json_local = batch_instance.read_input(pm5_json)
    clinvar_decisions_local = batch_instance.read_input(clinvar_decisions)

    # create a new folder, move the files into it
    # unpack all the already compressed HTs into it
    # the compress once, containing all files for distribution
    job.command(
        f"""
        mkdir clinvarbitration_data
        mv {pm5_json_local} clinvarbitration_data/clinvar_decisions.pm5.json
        tar -xf {pm5_ht_local} -C clinvarbitration_data
        tar -xf {clinvar_decisions_local} -C clinvarbitration_data
        tar -czf {job.output} \
            clinvarbitration_data/clinvar_decisions.pm5.json \
            clinvarbitration_data/clinvar_decisions.ht \
            clinvarbitration_data/clinvar_decisions.pm5.ht
    """,
    )
    batch_instance.write_output(job.output, output)

    return job
