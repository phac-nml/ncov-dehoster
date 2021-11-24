// Nanopore dehosting workflows

// Enable dsl2
nextflow.enable.dsl = 2

// Import modules
// From nanostripper pipeline
include {nanostripper} from '../modules/nanopore_nanostripper.nf' 
include {guppyBasecallerGPU} from '../modules/nanopore_nanostripper.nf'
include {guppyBasecallerCPU} from '../modules/nanopore_nanostripper.nf'
include {combineFastq} from '../modules/nanopore_nanostripper.nf'
include {fastqSizeSelection} from '../modules/nanopore_nanostripper.nf'
include {fastqDemultiplex} from '../modules/nanopore_nanostripper.nf'
include {combineFast5Barcodes} from '../modules/nanopore_nanostripper.nf'
include {regenerateFast5s} from '../modules/nanopore_nanostripper.nf'
include {generateSimpleSequencingSummary} from '../modules/nanopore_nanostripper.nf'
include {combineCSVs} from '../modules/nanopore_nanostripper.nf'

// From minimap2 pipeline
include {generateMinimap2Index} from '../modules/nanopore_minimap2.nf'
include {guppyplexSizeSelection} from '../modules/nanopore_minimap2.nf'
include {compositeMapping} from '../modules/nanopore_minimap2.nf'
include {removeHumanReads} from '../modules/nanopore_minimap2.nf'
include {generateFastqFiles} from '../modules/nanopore_minimap2.nf'


// Workflow Nanostripper
workflow nanoporeNanostripperDehosting {
    take:
      ch_fast5
      ch_HumanReference
      ch_CovidReference
    
    main:

      if ( params.human_minimap2_index ) {
        Channel.fromPath( "${params.human_minimap2_index}")
                        .set{ ch_human_index }

        nanostripper(ch_fast5
                    .combine(ch_CovidReference)
                    .combine(ch_human_index))

      } else {
        nanostripper(ch_fast5
                    .combine(ch_CovidReference)
                    .combine(ch_HumanReference))
      }

        nanostripper.out.dehostedFast5.collect()
                                      .set{ ch_dehosted_only_fast5 }

        combineFast5Barcodes(ch_dehosted_only_fast5)


      if ( params.guppyCPU ) {

        // Need guppy to de-multiplex
        // Save guppy gpu for basecalling as its the slowest step and we have limited number
        // Can make it so that it can be only GPU guppy but for the moment keep it like this
        if ( params.guppyGPU ) {
          guppyBasecallerGPU(ch_dehosted_only_fast5)

          guppyBasecallerGPU.out
                            .set{ ch_basecalled_fastqs }

        } else {
          guppyBasecallerCPU(nanostripper.out.dehostedFast5)

          guppyBasecallerCPU.out.collect()
                            .set{ ch_basecalled_fastqs }
        }

        // Back to the same processes after basecalling
        combineFastq(ch_basecalled_fastqs)

        fastqSizeSelection(combineFastq.out)

        fastqDemultiplex(fastqSizeSelection.out)

        regenerateFast5s(fastqDemultiplex.out.barcodes.flatten(),
                          combineFast5Barcodes.out)

        generateSimpleSequencingSummary(regenerateFast5s.out.fast5_pass.collect())

        combineCSVs(regenerateFast5s.out.csv.collect())

      } else {
        println('WARNING: dehosted fast5 files cannot be basecalled without a specified guppy environment, dehosting fast5 files only and then exiting')
      }

}

// Workflow Minimap2 fastq files
workflow nanoporeMinimap2Dehosting {
  take:
      ch_fastq
      ch_HumanReference
      ch_CovidReference
      directoryIn
    
    main:

    // If given a reference utilize it, otherwise make it
    if ( params.composite_minimap2_index ) {
      Channel.fromPath( "${params.composite_minimap2_index}")
                        .set{ ch_CompReference }
    } else {
      generateMinimap2Index(ch_HumanReference, 
                                ch_CovidReference)
      generateMinimap2Index.out
                           .set{ ch_CompReference }
    }
    
    guppyplexSizeSelection(ch_fastq)
    
    compositeMapping(guppyplexSizeSelection.out
                                           .combine(ch_CompReference))

    removeHumanReads(compositeMapping.out)

    generateFastqFiles(removeHumanReads.out)

    if ( directoryIn ) {
      if ( params.fast5_directory ) {
        println('Directory passed correctly')
      }
    }
}