// Nanopore dehosting workflow

// Enable dsl2
nextflow.preview.dsl = 2

// Import modules
include {nanostripper} from '../modules/nanopore.nf' 
include {guppyBasecaller} from '../modules/nanopore.nf'


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

      if ( params.guppy ) {
        guppyBasecaller()

      } else {
        println('WARNING: dehosted fast5 files cannot be basecalled without a specified guppy, dehosting fast5 files and then exiting')
      }

}
