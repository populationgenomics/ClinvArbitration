from typing import TYPE_CHECKING

from cpg_utils import Path, config, hail_batch

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def clinvarbitration_nextflow(
    output_root: str,
) -> 'BashJob':
    """
    runs the whole process using nextflow, in a single stage
    more a proof of concept than actually useful in the CPG deployment
    """

    batch_instance = hail_batch.get_batch('Run ClinvArbitration')

    # read in a reference genome
    ref_fa = batch_instance.read_input(config.config_retrieve(['workflow', 'ref_fa']))

    # make a new job
    job = batch_instance.new_bash_job('Run ClinvArbitration Nextflow')

    job.image(config.config_retrieve(['workflow', 'driver_image']))

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

    gff3 = batch_instance.read_input(config.config_retrieve(['references', 'ensembl_113', 'gff3']))

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
    batch_instance.write_output(job.output, output_root)
    return job
