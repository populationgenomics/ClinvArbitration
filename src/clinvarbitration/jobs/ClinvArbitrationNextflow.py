from typing import TYPE_CHECKING

from cpg_utils.config import config_retrieve
from cpg_utils.hail_batch import get_batch


if TYPE_CHECKING:
    from hailtop.batch.job import BashJob
    from cpg_utils import Path


def clinvarbitration_nextflow(
    output_root: 'Path',
) -> 'BashJob':
    """
    runs the whole process using nextflow, in a single stage
    more a proof of concept than actually useful
    Args:
        output_root ():

    Returns:

    """

    # read in a reference genome
    ref_fa = get_batch().read_input(config_retrieve(['workflow', 'ref_fa']))

    # make a new job
    job = get_batch().new_bash_job('Run ClinvArbitration Nextflow')

    job.image(config_retrieve(['workflow', 'driver_image']))

    # set some resource params
    job.storage('10Gi').memory('highmem').cpu(2)

    # set up one sprawling output group for all workflow results
    job.declare_resource_group(
        output={
            'submission_raw.txt.gz': '{root}/submission_summary.txt.gz',
            'variant_raw.txt.gz': '{root}/variant_summary.txt.gz',
            'clinvar_decisions.json': '{root}/clinvar_decisions.json',
            'ht.tar.gz': '{root}/clinvar_decisions.ht.tar.gz',
            'pm5.ht.tar.gz': '{root}/clinvar_decisions.pm5.ht.tar.gz',
            'pm5.json': '{root}/clinvar_decisions.pm5.json',
            'vcf.bgz': '{root}/clinvar_decisions.vcf.bgz',
            'vcf.bgz.tbi': '{root}/clinvar_decisions.vcf.bgz.tbi',
            'unfiltered.vcf.bgz': '{root}/clinvar_decisions.unfiltered.vcf.bgz',
            'unfiltered.vcf.bgz.tbi': '{root}/clinvar_decisions.unfiltered.vcf.bgz.tbi',
            'annotated.tsv': '{root}/clinvar_decisions.annotated.tsv',
            'release.tar.gz': '{root}/clinvar_decisions.release.tar.gz',
        },
    )

    gff3 = get_batch().read_input(config_retrieve(['references', 'ensembl_113', 'gff3']))

    # nextflow go brrrr
    job.command(
        f"""
        nextflow \
            -c nextflow/nextflow.config \
            run nextflow/clinvarbitration.nf \
            --ref_fa {ref_fa} \
            --output_dir {job.output} \
            --gff3 {gff3}
        """,
    )

    # copy the outputs back, in one smooooooth motion
    get_batch().write_output(job.output, output_root)
