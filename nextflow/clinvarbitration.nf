#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { AnnotateCsqWithBcftools } from './modules/clinvarbitration/AnnotateCsqWithBcftools/main'
include { LocaliseVariants } from './modules/clinvarbitration/LocaliseVariants/main'
include { LocaliseSubmissions } from './modules/clinvarbitration/LocaliseSubmissions/main'
include { MakePm5TableFromAnnotations } from './modules/clinvarbitration/MakePm5TableFromAnnotations/main'
include { PackageForRelease } from './modules/clinvarbitration/PackageForRelease/main'
include { ResummariseRawSubmissions } from './modules/clinvarbitration/ResummariseRawSubmissions/main'

params.publish_mode = 'copy'

check_params()

workflow {
    // localise both input files (in parallel, they're chunky)
    LocaliseVariants()
    LocaliseSubmissions()

    // reinterpret the results using altered heuristics
    ResummariseRawSubmissions(
        LocaliseVariants.out.variants,
        LocaliseSubmissions.out.submissions,
    )

    // get and check the reference genome. /Users/mwelland/binfx/refgenomes/Homo_sapiens_assembly38_masked.fasta
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
