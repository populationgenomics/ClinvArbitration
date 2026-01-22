#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { AnnotateCsqWithBcftools } from './modules/AnnotateCsqWithBcftools/main'
include { DownloadClinVarFiles } from './modules/DownloadClinVarFiles/main'
include { MakePm5TableFromAnnotations } from './modules/MakePm5TableFromAnnotations/main'
include { PackageForRelease } from './modules/PackageForRelease/main'
include { ResummariseRawSubmissions } from './modules/ResummariseRawSubmissions/main'

params.publish_mode = 'copy'

check_params()

workflow {
	// These files are all downloaded by the script `data/download_files.sh`, run that prior to starting the workflow
	ch_gff3 = channel.fromPath(params.gff3, checkIfExists:true)
    ch_ref_fa = channel.fromPath(params.ref_fa, checkIfExists: true)

    def current_month = new java.util.Date().format('yyyy-MM')
    String subfile = "${params.data}/submissions_${current_month}.txt.gz"
    String varfile = "${params.data}/variants_${current_month}.txt.gz"

    if (file(subfile).exists() && file(varfile).exists()) {
        ch_clinvar_sub = Channel.fromPath(subfile)
        ch_variants = Channel.fromPath(varfile)
    } else {
        println "Attempting to download ClinVar raw data, requires internet connection."
        println "If this step fails, try re-running download_files.sh in the `data` directory."

        // this step requires an internet connection, which may be problematic at some sites
        DownloadClinVarFiles()
        ch_clinvar_sub = DownloadClinVarFiles.out.submissions
        ch_variants = DownloadClinVarFiles.out.variants
    }

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
