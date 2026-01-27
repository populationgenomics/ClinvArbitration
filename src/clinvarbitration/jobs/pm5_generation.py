from typing import TYPE_CHECKING

from cpg_utils import hail_batch

from clinvarbitration.cpg_internal.utils import make_me_a_job

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def generate_pm5_data(
    annotated_snvs: str,
    output_folder: str,
) -> 'BashJob':
    """Generate PM5 data (index pathogenic missense variants by codon/transcript)."""

    batch_instance = hail_batch.get_batch('Run ClinvArbitration')

    annotated_snvs_local = batch_instance.read_input(annotated_snvs)

    job = make_me_a_job('Pm5TableGeneration').storage('10G')

    # write both HT and TSV outputs to the same root location
    job.command(f"""
        python3 -m clinvarbitration.scripts.clinvar_by_codon \\
            -i {annotated_snvs_local} \\
            -o ${{BATCH_TMPDIR}}/clinvar_decisions.pm5
    """)

    # compress the HT, move the TSV, and copy everything out in a single command
    job.command(f"""
        mv ${{BATCH_TMPDIR}}/clinvar_decisions.pm5.ht clinvar_decisions.pm5.ht
        tar -cf clinvar_decisions.pm5.ht.tar clinvar_decisions.pm5.ht
        gcloud storage cp ${{BATCH_TMPDIR}}/clinvar_decisions.pm5.tsv clinvar_decisions.pm5.ht.tar {output_folder}
        """,
    )

    return job
