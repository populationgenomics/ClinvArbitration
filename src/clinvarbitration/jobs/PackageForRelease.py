from typing import TYPE_CHECKING

from cpg_utils import Path, config, hail_batch

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def package_data_for_release(
    annotated_tsv: Path,
    pm5_json: Path,
    pm5_ht: Path,
    clinvar_decisions: Path,
    output: Path,
) -> 'BashJob':
    """
    Localise all the previously generated data into a folder
    tarball it, and write out as a single file

    Args:
        annotated_tsv (Pathlike): the annotated SNVs in TSV format
        pm5_json (Pathlike): the PM5 JSON file
        pm5_ht (Pathlike): the PM5 Hail Table
        clinvar_decisions (Pathlike): the ClinVar decisions Hail Table
        output (Pathlike): path to write output file to
    """
    batch_instance = hail_batch.get_batch('Run ClinvArbitration')

    job = batch_instance.new_bash_job('PackageForRelease')
    job.image(config.config_retrieve(['workflow', 'driver_image'])).storage('10G')

    annotated_tsv_local = batch_instance.read_input(annotated_tsv)
    pm5_ht_local = batch_instance.read_input(pm5_ht)
    pm5_json_local = batch_instance.read_input(pm5_json)
    clinvar_decisions_local = batch_instance.read_input(clinvar_decisions)

    # create a new folder, move the files into it
    # unpack all the already compressed HTs into it
    # the compress once, containing all files for distribution
    job.command(
        f"""
        mkdir clinvarbitration_data
        mv {annotated_tsv_local} clinvarbitration_data/clinvar_decisions.annotated.tsv
        mv {pm5_json_local} clinvarbitration_data/clinvar_decisions.pm5.json
        tar -xf {pm5_ht_local} -C clinvarbitration_data
        tar -xf {clinvar_decisions_local} -C clinvarbitration_data
        tar -czf {job.output} \
            clinvarbitration_data/clinvar_decisions.annotated.tsv \
            clinvarbitration_data/clinvar_decisions.pm5.json \
            clinvarbitration_data/clinvar_decisions.ht \
            clinvarbitration_data/clinvar_decisions.pm5.ht
    """,
    )
    batch_instance.write_output(job.output, output)

    return job
