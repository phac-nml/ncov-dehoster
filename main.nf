#!/usr/bin/env nextflow

// enable dsl2
nextflow.preview.dsl = 2

// Modules
include copyReference from './modules/minimap.nf'
include minimap2 from './modules/minimap.nf'
include samtoolsFlagstat from './modules/minimap.nf'
include removeMappedReads from './modules/minimap.nf'
include extractFastq from './modules/minimap.nf'


if ( !params.directory ) {
    println('Please supply a directory containing fastq files with directory.')
    System.exit(1)
}

// Main
workflow {

    Channel.fromPath( "${params.directory}/*.fastq", type: 'file', maxDepth: 1 )
                        .set{ ch_fastq }

    Channel.fromPath( "${params.human_ref}")
                        .set{ ch_href }

    Channel.fromPath( "${params.cov2019_ref}")
                        .set{ ch_cref }

    copyReference(ch_href)

    minimap2(ch_fastq
                .combine(copyReference.out)
                .combine(ch_cref))

    samtoolsFlagstat(minimap2.out)

    removeMappedReads(minimap2.out)

    extractFastq(removeMappedReads.out)
}