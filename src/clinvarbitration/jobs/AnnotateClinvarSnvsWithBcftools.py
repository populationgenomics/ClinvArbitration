from typing import TYPE_CHECKING

from cpg_utils import Path, config, hail_batch

if TYPE_CHECKING:
    from hailtop.batch.job import BashJob


def annotate_clinvar_snvs(
    snv_vcf: Path,
    output: Path,
) -> 'BashJob':
    """Annotate the SNVs, reduce dataset to Missense only, write as a TSV."""

    batch_instance = hail_batch.get_batch('Run ClinvArbitration')

    # localise the input file
    snv_vcf_local = batch_instance.read_input(snv_vcf)

    # need a genome reference - mandatory argument
    ref_fa = batch_instance.read_input(config.config_retrieve(['workflow', 'ref_fa']))

    # read in the gene model file
    gff3 = batch_instance.read_input(config.config_retrieve(['references', 'ensembl_113', 'gff3']))

    job = batch_instance.new_bash_job('AnnotateClinvarSnvsWithBcftools', attributes={'tool': 'bcftools'})
    job.image(config.config_retrieve(['workflow', 'driver_image']))

    # -g is the GFF3 file, -f is the reference fasta
    # --local-csq is required to apply non-phase aware annotation
    # --force is required to use annotations without phase data
    job.command(f'bcftools csq --force --local-csq -f {ref_fa} -g {gff3} {snv_vcf_local} -o out.vcf')

    # split the bcftools CSQ fields, filter to missense, and write out a tab-delimited file
    # -d - duplicate, writes each transcript onto a new line
    # -s :missense - only keep Consequence==missense variants
    # -f - format string - tab delimited, Transcript, Amino Acid Change, ClinVar allele ID, ClinVar gold stars
    job.command(
        f"""
        bcftools \\
            +split-vep out.vcf \\
            -d \\
            -s :missense \\
            -f "%transcript\t%amino_acid_change\t%allele_id\t%gold_stars\n" \\
            > {job.output}
        """
    )

    batch_instance.write_output(job.output, output)

    return job
