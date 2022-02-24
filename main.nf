#!/usr/bin/env nextflow

// enable dsl2
nextflow.enable.dsl = 2

// Modules
include {helpStatement} from './modules/help.nf'

// subworkflows
include {illuminaDehosting} from './workflows/illumina_dehosting.nf'
include {nanoporeNanostripperDehosting} from './workflows/nanopore_dehosting.nf'
include {nanoporeMinimap2Dehosting} from './workflows/nanopore_dehosting.nf'

// Print Help
if ( params.help ){
    helpStatement()
    exit 0
}

// Checking that everything is found before starting process
if ( params.illumina ) {
    if ( !params.directory ) {
        println("Please specify a directory containing paired fastq reads with --directory <path/to/fastqs>")
        System.exit(1)
    }
} else if ( params.nanopore ) {
    // Pick a pipeline and validate inputs for it
    if ( !params.nanostripper && !params.minimap2 ) {
        println("Please specify the tool to dehost nanopore data")
        System.exit(1)
    } else if ( params.nanostripper && params.minimap2 ) {
        println("Please specify one of `--nanostripper` or `--minimap2`")
        System.exit(1)
    // Nanostripper Validation
    } else if ( params.nanostripper ) {
        if ( !params.fast5_directory ) {
            println("Please specify a directory containing barcoded fast5 files with --fast5_directory <path/to/BARCODES/>")
            System.exit(1)
        }
    // Minimap2 Validation
    } else {
        if ( !params.fastq_directory ) {
            println("Please specify a directory containing fastq files or barcoded fastq directories with --fastq_directory <path/to/fastqs>")
            System.exit(1)
        }
    }

    // Run name for sorting out nanopore outputs in default output directory
    if ( !params.run_name) {
        println("Please specify a run name for seperating out nanopore runs")
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
        
        // Single-end or paired-end run - We don't support single end Illumina at this time
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
        // Nanostripper with Fast5s input
        if ( params.nanostripper ) {
            // First check if barcoded
            barcodedFast5 = file("${params.fast5_directory}/*{barcode,unclassified}*", type: 'dir', maxDepth: 1)
            nonBarcodedFast5 = file("${params.fast5_directory}/*.fast5", type: 'file', maxDepth: 1)

            // Use barcode to parallelize running if there are any
            // Doesn't like softlinked directories, checks that we have files
            // If not done, errors on blank barcoded directories
            if ( barcodedFast5 ) {
                Channel.fromPath( barcodedFast5 )
                    .filter{ d ->
                                def count = 0
                                for (x in d.listFiles()) {
                                    if (x.isFile()) {
                                        count += 1
                                    }
                                }
                                count > 0
                    }.set{ ch_fast5 }

                nanoporeNanostripperDehosting(ch_fast5, ch_HumanReference, ch_CovidReference)

            } else if ( nonBarcodedFast5 ) {
                println('Non Barcoded Fast5 files unavailable to dehost')
                System.exit(1)

            } else {
                println('Unable to figure out input')
                System.exit(1)
            }

        // Minimap2 with Fastqs input
        } else {
            // First check if barcoded or not, need to 
            barcodedFastq = file("${params.fastq_directory}/*{barcode,unclassified}*", type: 'dir', maxDepth: 1)
            nonBarcodedFastq = file("${params.fastq_directory}/*.fastq*", type: 'file', maxDepth: 1)
            if ( barcodedFastq ) {
                Channel.fromPath( barcodedFastq )
                    .filter{ d ->
                                def count = 0
                                for (x in d.listFiles()) {
                                    if (x.isFile()) {
                                        count += 1
                                    }
                                }
                                count > 0
                    }.set{ ch_fastq }
            
            } else if ( nonBarcodedFastq ) {
                Channel.fromPath( nonBarcodedFastq )
                    .set{ ch_fastq }
            } else {
                println('Unable to figure out input')
                System.exit(1)
            }
            nanoporeMinimap2Dehosting(ch_fastq, ch_HumanReference, ch_CovidReference)
        }
    }
}
