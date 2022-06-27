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

// Workflow
workflow illuminaDehosting {
    take:
      ch_fastqs
      ch_HumanReference
      ch_CovidReference

    main:

    // Always need to make the composite reference, even if an index is given
    // This might lead to issues if the composite index given does match the generated reference but we will watch for that and re-visit it later
    generateCompositeReference(ch_HumanReference, 
                                ch_CovidReference)

    if ( params.composite_bwa_index ){
      grabCompositeIndex("${params.composite_bwa_index}")
      grabCompositeIndex.out
              .set{ ch_index }
    } else {
      indexCompositeReference(generateCompositeReference.out)
      indexCompositeReference.out.collect()
              .set{ ch_index }
    }

    compositeMappingBWA(ch_fastqs
                        .combine(generateCompositeReference.out),
                      ch_index)

    dehostBamFiles(compositeMappingBWA.out.bam)

    generateDehostedReads(dehostBamFiles.out.bam)

    combineCSVs(dehostBamFiles.out.csv.collect())
}
