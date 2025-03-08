#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { AnnotateCsqWithBcftools } from './modules/clinvarbitration/AnnotateCsqWithBcftools/main'
include { LocaliseVariants } from './modules/clinvarbitration/LocaliseVariants/main'
include { LocaliseSubmissions } from './modules/clinvarbitration/LocaliseSubmissions/main'
include { MakePm5TableFromAnnotations } from './modules/clinvarbitration/MakePm5TableFromAnnotations/main'
include { PackageForRelease } from './modules/clinvarbitration/PackageForRelease/main'
include { ResummariseRawSubmissions } from './modules/clinvarbitration/ResummariseRawSubmissions/main'

// since we're containerising these need to be duplicates, not symlinks, for when the container closes
params.publish_mode = 'copy'

check_params()

workflow {
    // localise both input files (in parallel, they're chunky) - check if they already exist
    // this just lets us resume the main processing steps instead of re-downloading this upon failure
    if (file(params.sub_sum_output).exists()) {
        println 'Re-using the localised SubmissionSummary.txt.gz file!'
        ch_sub_summary = channel.fromPath(params.sub_sum_output)
    }
    else {
        LocaliseSubmissions()
        ch_sub_summary = LocaliseSubmissions.out.submissions
    }
    if (file(params.var_sum_output).exists()) {
        println 'Re-using the localised VariantSummary.txt.gz file!'
        ch_var_summary = channel.fromPath(params.var_sum_output)
    }
    else {
        LocaliseVariants()
        ch_var_summary = LocaliseVariants.out.variants
    }

    // reinterpret the results using altered heuristics
    ResummariseRawSubmissions(
        ch_var_summary,
        ch_sub_summary,
    )

    // get and check the reference genome
    ch_ref_fa = channel.fromPath(params.ref_fa, checkIfExists: true)

    // annotate the SNV VCF using BCFtools
    AnnotateCsqWithBcftools(
        ResummariseRawSubmissions.out.vcf,
        ch_ref_fa
    )

    MakePm5TableFromAnnotations(
        AnnotateCsqWithBcftools.out.tsv
    )

    PackageForRelease(
        ResummariseRawSubmissions.out.ht_tar,
        MakePm5TableFromAnnotations.out.ht_tar,
        MakePm5TableFromAnnotations.out.json,
    )
}


// this triggers if "--help" is used
def check_params() {
    if( params.remove('help') ) {
        params.each{ k, v -> println "params.${k.padRight(25)} = ${v}" }
        exit 0
    }
}
