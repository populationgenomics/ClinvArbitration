process PackageForRelease {

    publishDir params.output_dir, mode: 'copy'

    input:
        path decisions_vcf
        path decisions_vcf_idx
        path decisions_tsv
        path pm5_json
        path pm5_tsv

    output:
        path "clinvar_decisions.release.tar.gz"

    script:
        // create a new folder, and compress everything together
        """
        mkdir clinvarbitration_data
        cp "${decisions_vcf}" "${decisions_vcf_idx}" "${decisions_tsv}" "${pm5_json}" "${pm5_tsv}" clinvarbitration_data/
        tar -czf clinvar_decisions.release.tar.gz clinvarbitration_data
        """
}
