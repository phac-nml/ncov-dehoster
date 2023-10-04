// Illumina dehosting workflow
// Import modules
include {
  generateCompositeReference;
  grabCompositeIndex;
  indexCompositeReference;
  compositeMappingBWA;
  dehostBamFiles;
  generateDehostedReads;
  combineCSVs
} from '../modules/illumina.nf'

include {
  seqtkRandomSubsample;
  samtoolsAmpliconDownsample
} from '../modules/general.nf'
include { outputVersions } from '../modules/versions.nf'

// Workflow
workflow illuminaDehosting {
    take:
      ch_fastqs
      ch_HumanReference
      ch_CovidReference

    main:

    // Setup tool version tracking - based on NF-Core's process
    ch_versions = Channel.empty()

    // Always need to make the composite reference, even if an index is given
    // This might lead to issues if the composite index given does match the generated reference but we will watch for that and re-visit it later
    generateCompositeReference(ch_HumanReference, 
                                ch_CovidReference)

    if ( params.composite_bwa_index ) {
      grabCompositeIndex("${params.composite_bwa_index}")
      grabCompositeIndex.out
              .set{ ch_index }
    } else {
      indexCompositeReference(generateCompositeReference.out)
      indexCompositeReference.out.index.collect()
              .set{ ch_index }

      ch_versions = ch_versions.mix(indexCompositeReference.out.versions)
    }

    compositeMappingBWA(ch_fastqs
                        .combine(generateCompositeReference.out),
                      ch_index)
    ch_versions = ch_versions.mix(compositeMappingBWA.out.versions.first())

    dehostBamFiles(compositeMappingBWA.out.bam)
    ch_versions = ch_versions.mix(dehostBamFiles.out.versions.first())

    // Downsampling BAM with Amplicons
    if ( params.downsample && params.downsample_amplicons ) {
      samtoolsAmpliconDownsample(
        dehostBamFiles.out.bam,
        file(params.downsample_amplicons, checkIfExists: true),
        'illumina',
        params.downsample_count,
        params.downsample_seed
      )

      ch_versions = ch_versions.mix(samtoolsAmpliconDownsample.out.versions.first())
      ch_dehosted_bam = samtoolsAmpliconDownsample.out.bam
    } else {
      ch_dehosted_bam = dehostBamFiles.out.bam
    }

    generateDehostedReads(ch_dehosted_bam)
    ch_versions = ch_versions.mix(generateDehostedReads.out.versions.first())

    // Downsampling Final Fastqs
    if ( params.downsample && !params.downsample_amplicons ) {
      seqtkRandomSubsample(generateDehostedReads.out.dehosted_fastq, params.downsample_count)
      ch_versions = ch_versions.mix(seqtkRandomSubsample.out.versions.first())
    }

    combineCSVs(dehostBamFiles.out.csv.collect())
    ch_versions = ch_versions.mix(combineCSVs.out.versions)

    // Version Tracking Output
    outputVersions(ch_versions.collectFile(name: 'tool_versions.yml'))
}
