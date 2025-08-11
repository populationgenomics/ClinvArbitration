#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { AnnotateCsqWithBcftools } from './modules/clinvarbitration/AnnotateCsqWithBcftools/main'
include { MakePm5TableFromAnnotations } from './modules/clinvarbitration/MakePm5TableFromAnnotations/main'
include { PackageForRelease } from './modules/clinvarbitration/PackageForRelease/main'
include { ResummariseRawSubmissions } from './modules/clinvarbitration/ResummariseRawSubmissions/main'

params.publish_mode = 'copy'

check_params()

workflow {
	// These files are all downloaded by the script `data/download_files.sh`, run that prior to starting the workflow
	ch_gff3 = channel.fromPath(params.gff3, checkIfExists:true)
	ch_submissions = channel.fromPath(params.clinvar_submissions, checkIfExists: true)
	ch_variants = channel.fromPath(params.clinvar_variants, checkIfExists: true)
    ch_ref_fa = channel.fromPath(params.ref_fa, checkIfExists: true)

    // reinterpret the results using altered heuristics
    ResummariseRawSubmissions(
        ch_variants,
        ch_submissions,
    )

    // annotate the SNV VCF using BCFtools
    AnnotateCsqWithBcftools(
        ResummariseRawSubmissions.out.vcf,
        ch_ref_fa,
        ch_gff3,
    )

    MakePm5TableFromAnnotations(
        AnnotateCsqWithBcftools.out
    )

    PackageForRelease(
        ResummariseRawSubmissions.out.ht,
        ResummariseRawSubmissions.out.tsv,
        MakePm5TableFromAnnotations.out.ht,
        MakePm5TableFromAnnotations.out.tsv,
    )
}


// this triggers if "--help" is used
def check_params() {
    if( params.remove('help') ) {
        params.each{ k, v -> println "params.${k.padRight(25)} = ${v}" }
        exit 0
    }
}
