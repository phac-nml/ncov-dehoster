// Nanopore dehosting workflow

// Enable dsl2
nextflow.preview.dsl = 2

// Import modules
include {nanostripper} from '../modules/nanopore.nf' 
include {guppyBasecallerGPU} from '../modules/nanopore.nf'
include {guppyBasecallerCPU} from '../modules/nanopore.nf'
include {combineFastq} from '../modules/nanopore.nf'
include {fastqSizeSelection} from '../modules/nanopore.nf'
include {fastqDemultiplex} from '../modules/nanopore.nf'
include {combineFast5Barcodes} from '../modules/nanopore.nf'
include {regenerateFast5s} from '../modules/nanopore.nf'


// Workflow
workflow nanoporeDehosting {
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

      } else {
        println('WARNING: dehosted fast5 files cannot be basecalled without a specified guppy environment, dehosting fast5 files and then exiting')
      }

}
