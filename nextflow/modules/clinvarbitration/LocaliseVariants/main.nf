
process LocaliseVariants {

    publishDir params.output_dir, mode: 'copy'

    output:
        path "variant_summary.txt.gz", emit: "variants"

    """
    wget -O "variant_summary.txt.gz" "${params.variants}"
    """
}
