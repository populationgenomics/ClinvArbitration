
process MakePm5TableFromAnnotations {

    publishDir params.output_dir, mode: 'copy'

    input:
        path annotated_snv

    output:
        path "clinvar_decisions.pm5.ht", emit: "ht"
        path "clinvar_decisions.pm5.tsv", emit: "tsv"

    """
    python3 -m clinvarbitration.scripts.clinvar_by_codon \
        -i "${annotated_snv}" \
        -o clinvar_decisions.pm5
    """
}
