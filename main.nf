#!/usr/bin/env nextflow

// enable dsl2
nextflow.preview.dsl = 2

// Modules
// None at the moment, add some here if need to pre-process reads (like find no inputs etc)

// subworkflows
include {illuminaDehosting} from './workflows/illumina_dehosting.nf'
//include {nanoporeDehosting} from './workflows/nanopore_dehosting.nf'


if ( !params.directory ) {
    println('Please supply a directory containing fastq files with directory.')
    System.exit(1)
}


// Checking everything is found
if ( params.illumina ) {
    if ( !params.directory ) {
        println("Please specify a directory containing single-end or paired fastq files with --directory <path/to/fastqs>")
        System.exit(1)
    }
} else if ( params.nanopore ) {
    if ( !params.directory ) {
        println("Please specify a directory containing fast5 files with --directory <path/to/fast5s>")
        System.exit(1)
    }
} else {
    println('Please specify either --illumina or --nanopore for the type of sequencing data that will be dehosted')
    System.exit(1)
}

// Main
workflow {

    Channel.fromPath( "${params.human_ref}")
                        .set{ ch_HumanReference }

    Channel.fromPath( "${params.cov2019_ref}")
                        .set{ ch_CovidReference }

    if ( params.illumina ) {
        
        // Single-end or paired-end run
        if ( params.single_end ){
            Channel.fromPath( "${params.directory}/*.fastq", type: 'file', maxDepth: 1 )
                        .set{ ch_fastqs }
            
            println('Single end reads pipeline is not yet available')
            System.exit(1)

        } else {
            Channel.fromFilePairs( params.fastqpaths, flat: true, maxDepth: 1)
                        .set{ ch_fastqs }

            illuminaDehosting(ch_fastqs, ch_HumanReference, ch_CovidReference)
        }
    
    } else if ( params.nanopore ) {
        println('Nanopore pipeline not ready yet. This feature will be available soon (hopefully)')
        System.exit(1)

    } else {
        println('Please specify either --illumina or --nanopore for the type of data that will be dehosted')
        System.exit(1)
    }

}
