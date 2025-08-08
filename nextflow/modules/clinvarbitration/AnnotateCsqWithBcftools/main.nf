
process AnnotateCsqWithBcftools {

    publishDir params.output_dir

    input:
        // the two input files from ClinVar
        path vcf
        path ref_fa
        path gff3

    output:
        path "clinvar_decisions.annotated.tsv"

    """
    bcftools csq \
        -f "${ref_fa}" \
        -g "${gff3}" \
        "${vcf}" |
    bcftools +split-vep \
        -d \
        -s :missense \
        -f "%transcript\t%amino_acid_change\t%allele_id\t%gold_stars\n" \
        - \
        > clinvar_decisions.annotated.tsv
    """
}
