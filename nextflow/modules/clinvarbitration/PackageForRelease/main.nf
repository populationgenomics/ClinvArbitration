
process PackageForRelease {

    publishDir params.output_dir, mode: 'copy'

    input:
        path decisions_ht_tar
        path pm5_ht_tar
        path pm5_json

    output:
        path "clinvar_decisions.release.tar.gz"

    // create a new folder, decompress the previous archives, and recompress everything together
    """
    mkdir clinvarbitration_data
    tar -xzf "${pm5_ht_tar}" -C clinvarbitration_data
    tar -xzf "${decisions_ht_tar}" -C clinvarbitration_data
    tar -czf clinvar_decisions.release.tar.gz \
        clinvarbitration_data/clinvar_decisions.ht \
        clinvarbitration_data/clinvar_decisions.pm5.ht \
        "${pm5_json}"
    """
}
