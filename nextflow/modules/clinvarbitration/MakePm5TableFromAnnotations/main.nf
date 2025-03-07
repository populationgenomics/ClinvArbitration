
process MakePm5TableFromAnnotations {

    publishDir params.output_dir, mode: 'copy'

    input:
        path annotated_snv

    output:
        path "clinvar_pm5.json", emit: "json"
        path "clinvar_pm5.ht.tar.gz", emit: "ht_tar"

    """
    pm5_table \
        -i "${annotated_snv}" -o clinvar_pm5

    tar -czf "clinvar_pm5.ht.tar.gz" clinvar_pm5.ht
    """
}
