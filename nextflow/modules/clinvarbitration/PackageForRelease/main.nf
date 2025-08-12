process PackageForRelease {

    publishDir params.output_dir, mode: 'copy'

    input:
        path decisions_ht
        path decisions_tsv
        path pm5_ht
        path pm5_tsv

    output:
        path "clinvar_decisions.release.tar.gz"

    // create a new folder, decompress the previous archives, and recompress everything together
    """
    mkdir clinvarbitration_data
    cp -r "${pm5_ht}" "${pm5_tsv}" "${decisions_ht}" "${decisions_tsv}" clinvarbitration_data/
    tar -czf clinvar_decisions.release.tar.gz clinvarbitration_data
    """
}
