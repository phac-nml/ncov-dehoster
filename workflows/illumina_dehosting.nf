// Illumina dehosting workflow

// Enable dsl2
nextflow.preview.dsl = 2

// Import modules
include {captureViralReads} from '../modules/illumina.nf'
include {captureHumanReads} from '../modules/illumina.nf'

include {copyReference} from '../modules/base_functions.nf'
include {extractFastq} from '../modules/base_functions.nf'

// Workflow
workflow illuminaDehosting {
    take:
      ch_fastqs
      ch_HumanReference
      ch_CovidReference

    main:

    copyReference(ch_HumanReference)

    captureViralReads(ch_fastqs
                        .combine(ch_CovidReference))

    captureHumanReads(ch_fastqs
                        .combine(copyReference.out))


}