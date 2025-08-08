process ResummariseRawSubmissions {

    publishDir params.output_dir, mode: 'copy'

    input:
        // the two input files from ClinVar
        path variant_summary
        path submission_summary

    output:
        path "clinvar_decisions.vcf.bgz", emit: "vcf"
        path "clinvar_decisions.vcf.bgz.tbi", emit: "vcf_idx"
        path "clinvar_decisions.ht", emit: "ht"

    // Generates
    // clinvar_decisions.vcf.bgz + index - VCF containing only pathogenic SNV entries, feeds into annotation
    // clinvar_decisions.ht - a Hail Table containing the summarised data entries
    """
    python3 -m clinvarbitration.scripts.resummarise_clinvar \
        -v "${variant_summary}" \
        -s "${submission_summary}" \
        -o "clinvar_decisions" \
        --assembly "${params.assembly}" \
        --minimal
    rm clinvar_decisions.json
    """
}
