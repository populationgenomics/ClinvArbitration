from typing import TYPE_CHECKING

from cpg_utils import config, hail_batch

from clinvarbitration.cpg_internal.utils import make_me_a_job

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def generate_new_summary(
    var_file: str,
    sub_file: str,
    output_root: str,
) -> 'BashJob':
    """Using the submission and variants data files, generate revised variant summaries."""
    batch_instance = hail_batch.get_batch()

    job = make_me_a_job('GenerateNewClinvarSummary').memory('highmem').cpu('2')

    if sites_to_blacklist := config.config_retrieve(['workflow', 'site_blacklist'], []):
        blacklist_sites = ' '.join(f'"{site}"' for site in sites_to_blacklist)
        blacklist_string = f' -b {blacklist_sites}'
    else:
        blacklist_string = ''

    var_file_local = batch_instance.read_input(var_file)
    sub_file_local = batch_instance.read_input(sub_file)

    job.command(f"""
        python3 -m clinvarbitration.scripts.resummarise_clinvar \\
        -v {var_file_local} \\
        -s {sub_file_local} \\
        {blacklist_string} -o ${{BATCH_TMPDIR}}/clinvar_decisions
    """)

    # don't tar from current location, we'll catch all the tmp pathing
    job.command(f"""
        mv ${{BATCH_TMPDIR}}/clinvar_decisions.ht clinvar_decisions.ht
        tar -cf clinvar_decisions.ht.tar clinvar_decisions.ht
        gcloud storage cp \\
            clinvar_decisions.ht.tar \\
            "${{BATCH_TMPDIR}}/clinvar_decisions.vcf.bgz*" \\
            ${{BATCH_TMPDIR}}/clinvar_decisions.tsv \\
            {output_root}
    """)

    return job
