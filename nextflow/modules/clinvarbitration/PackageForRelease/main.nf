
process PackageForRelease {

    publishDir params.output_dir, mode: 'copy'

    input:
        path decisions_ht_tar
        path pm5_ht_tar
        path pm5_json

    output:
        path "clinvarbitration_data.tar.gz"

    //
    """
    mkdir clinvarbitration_data
    tar -xzf "${pm5_ht_tar}" -C clinvarbitration_data
    tar -xzf "${decisions_ht_tar}" -C clinvarbitration_data
    tar -czf clinvarbitration_data.tar.gz \
        clinvarbitration_data/clinvar_decisions.ht \
        clinvarbitration_data/clinvar_pm5.ht \
        "${pm5_json}"
    """
}
