/// Nanopore Dehosting Workflows
// Import modules //
// From minimap2 pipeline
include {
  generateMinimap2Index;
  fastqSizeSelection_MM2;
  compositeMappingMM2;
  removeHumanReads;
  samtoolsSort ;
  regenerateFastqFiles;
  regenerateFast5s_MM2
} from '../modules/nanopore_minimap2.nf'

include {
  seqtkRandomSubsample ;
  samtoolsAmpliconDownsample ;
} from '../modules/general.nf'
include { outputVersions } from '../modules/versions.nf'

// From nanostripper pipeline
include {
  nanostripper;
  guppyBasecallerGPU;
  guppyBasecallerCPU;
  combineFastq;
  fastqSizeSelection_NS;
  fastqDemultiplex;
  combineFast5Barcodes;
  regenerateFast5s_NS;
  generateSimpleSequencingSummary;
  combineCSVs
} from '../modules/nanopore_nanostripper.nf'

// Workflow Minimap2 fastq files //
workflow nanoporeMinimap2Dehosting {
  take:
      ch_fastq
      ch_HumanReference
      ch_CovidReference
    
    main:
    // Setup tool version tracking - based on NF-Core's process
    ch_versions = Channel.empty()

    // If given a reference utilize it, otherwise make it
    if ( params.composite_minimap2_index ) {
      Channel.fromPath("${params.composite_minimap2_index}")
                        .set{ ch_CompReference }
    } else {
      generateMinimap2Index(ch_HumanReference, 
                                ch_CovidReference)
      generateMinimap2Index.out.index
                           .set{ ch_CompReference }

      ch_versions = ch_versions.mix(generateMinimap2Index.out.versions)
    }
    
    fastqSizeSelection_MM2(ch_fastq)
    ch_versions = ch_versions.mix(fastqSizeSelection_MM2.out.versions.first())
    
    compositeMappingMM2(fastqSizeSelection_MM2.out.fastq
                                           .combine(ch_CompReference))
    ch_versions = ch_versions.mix(compositeMappingMM2.out.versions.first())

    removeHumanReads(compositeMappingMM2.out.comp_bam)
    ch_versions = ch_versions.mix(removeHumanReads.out.versions.first())
    ch_dehosted_bam = removeHumanReads.out.bam

    // Downsampling BAM with Amplicons
    if ( params.downsample && params.downsample_amplicons ) {
      samtoolsAmpliconDownsample(
        removeHumanReads.out.bam,
        file(params.downsample_amplicons, checkIfExists: true),
        'nanopore',
        params.downsample_count,
        params.downsample_seed
      )

      ch_versions = ch_versions.mix(samtoolsAmpliconDownsample.out.versions.first())
      ch_dehosted_bam = samtoolsAmpliconDownsample.out.bam
    }

    samtoolsSort(ch_dehosted_bam)

    // Output either flat fastq directory or normal nanopore formatted output based on CL --flat arg
    regenerateFastqFiles(samtoolsSort.out.bam)
    regenerateFastqFiles.out.dehosted_fastq
                        .filter{ it[1].countFastq() >= params.min_read_count }
                        .set { ch_host_rm_fastq }

    ch_versions = ch_versions.mix(regenerateFastqFiles.out.versions.first())

    // Downsampling Fastq
    if ( params.downsample && !params.downsample_amplicons ) {
      seqtkRandomSubsample(ch_host_rm_fastq, params.downsample_count)

      ch_versions = ch_versions.mix(seqtkRandomSubsample.out.versions.first())
      ch_final_fastq = seqtkRandomSubsample.out.reads
    } else {
      ch_final_fastq = ch_host_rm_fastq
    }

    // If a fast5 directory is given, we can use the fastq files to regenerate dehosted fast5 files
    // This process is slow without a lot of computational support behind it however; so its optional
    if ( params.fast5_directory ) {
      Channel.fromPath( "${params.fast5_directory}")
                        .set{ ch_Fast5 }
      regenerateFast5s_MM2(ch_final_fastq.combine(ch_Fast5))

      ch_versions = ch_versions.mix(regenerateFast5s_MM2.out.versions.first())

      generateSimpleSequencingSummary(regenerateFast5s_MM2.out.dehosted_fast5.collect())
    }

    // Finally make CSV output
    combineCSVs(removeHumanReads.out.csv.collect())
    ch_versions = ch_versions.mix(combineCSVs.out.versions)

    // Version Tracking and Output
    outputVersions(ch_versions.collectFile(name: 'tool_versions.yml'))
}

// Workflow Nanostripper - !!!NOT MAINTAINED AT THE MOMENT!!!//
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
        // Guppy basecalling with CPU is insanely slow as a heads up so this will could take >12h
        } else {
          guppyBasecallerCPU(nanostripper.out.dehostedFast5)

          guppyBasecallerCPU.out.collect()
                            .set{ ch_basecalled_fastqs }
        }

        // Back to the same processes after basecalling
        combineFastq(ch_basecalled_fastqs)

        fastqSizeSelection_NS(combineFastq.out)

        fastqDemultiplex(fastqSizeSelection_NS.out)

        regenerateFast5s_NS(fastqDemultiplex.out.barcodes.flatten(),
                          combineFast5Barcodes.out)

        generateSimpleSequencingSummary(regenerateFast5s_NS.out.fast5_pass.collect())

        combineCSVs(regenerateFast5s_NS.out.csv.collect())

      } else {
        println('WARNING: dehosted fast5 files cannot be basecalled without a specified guppy environment, dehosting fast5 files only and then exiting')
      }
}
