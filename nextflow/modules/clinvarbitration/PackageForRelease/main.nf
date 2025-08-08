process PackageForRelease {

    publishDir params.output_dir, mode: 'copy'

    input:
        path decisions_ht
        path pm5_ht

    output:
        path "clinvar_decisions.release.tar"

    // create a new folder, decompress the previous archives, and recompress everything together
    """
    mkdir clinvarbitration_data
    mv "${pm5_ht}" clinvarbitration_data/
    mv "${decisions_ht}" clinvarbitration_data/
    tar -cf clinvar_decisions.release.tar clinvarbitration_data
    """
}
