
process AnnotateCsqWithBcftools {

    publishDir params.output_dir, mode: 'copy'

    input:
        // the two input files from ClinVar
        path vcf
        path ref_fa

    output:
        path "clinvar_decisions.annotated.tsv", emit: "tsv"

    """
    bcftools csq \
        -f "${ref_fa}" \
        -g "${params.gff3}" \
        "${vcf}" \
        -o temp_annotated_output.vcf
    bcftools +split-vep \
        temp_annotated_output.vcf \
        -d \
        -s :missense \
        -f "%transcript\t%amino_acid_change\t%allele_id\t%gold_stars\n" \
        > clinvar_decisions.annotated.tsv
    """
}
