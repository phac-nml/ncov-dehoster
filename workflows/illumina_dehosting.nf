// Illumina dehosting workflow

// Enable dsl2
nextflow.preview.dsl = 2

// Import modules
include {indexHumanReference} from '../modules/illumina.nf'
include {indexViralReference} from '../modules/illumina.nf'
include {captureViralReads} from '../modules/illumina.nf'
include {captureHumanReads} from '../modules/illumina.nf'
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

    Channel.fromPath("${params.human_bwa_index}")
              .set{ ch_index }

    indexHumanReference(ch_HumanReference, ch_index)
    indexViralReference(ch_CovidReference)

    captureViralReads(ch_fastqs
                        .combine(ch_CovidReference),
                      indexViralReference.out.collect())

    captureHumanReads(ch_fastqs
                        .combine(ch_HumanReference),
                      indexHumanReference.out.collect())

    dehostSamFiles(captureViralReads.out.sam
                    .join(captureHumanReads.out.sam, by: 0))

    generateDehostedReads(dehostSamFiles.out.sam)

    combineCSVs(dehostSamFiles.out.csv.collect())

}