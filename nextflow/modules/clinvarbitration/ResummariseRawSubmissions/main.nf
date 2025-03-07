
process ResummariseRawSubmissions {

    publishDir params.output_dir, mode: 'copy'

    input:
        // the two input files from ClinVar
        path variant_summary
        path submission_summary

    output:
        path "clinvar_decisions.json", emit: "json"
        path "clinvar_decisions.vcf.bgz", emit: "vcf"
        path "clinvar_decisions_unfiltered.vcf.bgz", emit: "vcf_all"
        path "clinvar_decisions.ht.tar.gz", emit: "ht_tar"

    // Generates
    // clinvar_decisions.json - a JSON file containing the summarised data entries, one json object per line
    // clinvar_decisions_unfiltered.vcf.bgz - a bgzipped VCF which can be used in VEP annotation (all decisions)
    // clinvar_decisions.vcf.bgz - a bgzipped VCF containing only pathogenic SNV entries
    // clinvar_decisions.ht.tar.gz - a Hail Table containing the summarised data entries, compressed
    """
    resummary \
        -v "${variant_summary}" \
        -s "${submission_summary}" \
        -o "clinvar_decisions" \
        --assembly "${params.assembly}" \
        --minimal
    tar -czf clinvar_decisions.ht.tar.gz clinvar_decisions.ht
    """
}
