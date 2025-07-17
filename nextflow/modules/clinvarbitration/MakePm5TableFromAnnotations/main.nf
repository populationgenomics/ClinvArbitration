
process MakePm5TableFromAnnotations {

    publishDir params.output_dir, mode: 'copy'

    input:
        path annotated_snv

    output:
        path "clinvar_decisions.pm5.json", emit: "json"
        path "clinvar_decisions.pm5.ht.tar.gz", emit: "ht_tar"

    """
    python3 -m clinvarbitration.scripts.clinvar_by_codon \
        -i "${annotated_snv}" \
        -o clinvar_decisions.pm5

    tar -czf "clinvar_decisions.pm5.ht.tar.gz" clinvar_decisions.pm5.ht
    """
}
