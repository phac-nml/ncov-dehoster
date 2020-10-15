// Illumina dehosting workflow

// Enable dsl2
nextflow.preview.dsl = 2

// Import modules
include {generateCompositeReference} from '../modules/illumina.nf'
include {grabCompositeIndex} from '../modules/illumina.nf'
include {indexCompositeReference} from '../modules/illumina.nf'
include {mapToCompositeIndex} from '../modules/illumina.nf'
include {dehostSamFiles} from '../modules/illumina.nf'
include {generateDehostedReads} from '../modules/illumina.nf'
include {combineCSVs} from '../modules/illumina.nf'

// Workflow
workflow illuminaDehosting {
    take:
      ch_fastqs
      ch_HumanReference
      ch_CovidReference

    main:

    generateCompositeReference(ch_HumanReference, 
                                ch_CovidReference)

    if ( params.composite_bwa_index ){
      grabCompositeIndex("${params.composite_bwa_index}")

      grabCompositeIndex.out
              .set{ ch_index }
    } else {
      indexCompositeReference(generateCompositeReference.out)

      indexCompositeReference.out
              .set{ ch_index }
    }

    mapToCompositeIndex(ch_fastqs
                        .combine(generateCompositeReference.out),
                      ch_index)

    dehostSamFiles(mapToCompositeIndex.out.sam)

    generateDehostedReads(dehostSamFiles.out.sam)

    combineCSVs(dehostSamFiles.out.csv.collect())
}
