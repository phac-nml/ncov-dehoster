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

include { seqtkSubsample } from '../modules/general.nf'
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

    generateDehostedReads(dehostBamFiles.out.bam)
    ch_versions = ch_versions.mix(generateDehostedReads.out.versions.first())

    // Downsampling
    if ( params.downsample ) {
      seqtkSubsample(generateDehostedReads.out.dehosted_fastq, params.downsample_count)
      ch_versions = ch_versions.mix(seqtkSubsample.out.versions.first())
    }

    combineCSVs(dehostBamFiles.out.csv.collect())
    ch_versions = ch_versions.mix(combineCSVs.out.versions)

    // Version Tracking Output
    outputVersions(ch_versions.collectFile(name: 'tool_versions.yml'))
}
