
process LocaliseSubmissions {

    publishDir params.output_dir, mode: 'copy'

    output:
        path "submission_summary.txt.gz", emit: "submissions"

    """
    wget -O "submission_summary.txt.gz" "${params.submissions}"
    """
}
