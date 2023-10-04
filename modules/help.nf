def helpStatement() {
  log.info"""
    Usage:
        nextflow run phac-nml/ncov-dehoster --human_ref [Path/to/Fasta] ( --illumina | --nanopore ( --nanostripper | --minimap2 )) [workflow-options]

    Overall Mandatory Arguments:
        --human_ref [file]              Path to human reference genome (ex. hg38.fa)

    Overall Optional Arguments:
        --cache [path]                  Path to cache directory containing conda environments
        --outdir [str]                  Name for the output directory (Default: 'results')
        --downsample                    Turn on random downsampling of reads with seqtk
        --downsample_count [i]          Integer approximate number of reads to downsample to (Default: 100,000)
        --downsample_seed [i]           Integer downsample seed to utilize (Default: 101)
        --downsample_amplicons [file]   Bed file to use to downsample amplicon regions with samtools instead of randomly with seqtk

    Profiles:
        -profile [conda,nml]            Configuration profile to use. Can use both if separated by a comma
            conda                         Utilize conda to control tool and dependency installation (uses mamba for environment install)
            nml                           NML specific profile to take advantage of NML cluster resources

    Pick a Pipeline:
        --illumina                      Run Illumina BWA_MEM based competitive host removal pipeline on paired fastq files
        --nanopore                      Set Nanopore data as input, MUST pick a sub pipeline
            --minimap2                    Run Nanopore host removal based on input fastq files and minimap2
            --nanostripper                Run Nanopore host removal based on input fast5 files, nanostripper, and guppy (NOT CURRENTLY MAINTAINED)

    Illumina Pipeline Arguments:
        Mandatory:
            --directory [path]          Path to directory containing paired fastq files to dehost

        Optional:
            --composite_bwa_index [path]    Path to directory containing BWA indexed composite genome (Note: still need human reference fasta)
            --keep_ref_id [str]             String name of reference contig ID to keep during host removal (Default: 'MN908947.3')
            --keep_min_map_quality [i]      Integer minimum mapping quality of COVID-19 reads to keep (Default: 60)
            --remove_min_map_quality [i]    Integer minimum mapping quality of HUMAN reads to remove (Default: 0)

        Example Command:
            nextflow run phac-nml/ncov-dehoster --human_ref ./hg38.fa -profile conda --illumina --directory ./paired_fastqs

                                        ------------------------------------------

    Minimap2 Pipeline Arguments:
        Mandatory:
            --fastq_directory [path]            Path to directory containing either barcoded directories OR individual fastq files
            --run_name [str]                    Seperate results based on a run name input

        Optional:
            --fast5_directory [path]            Path to fast5 directories associated with the data
            --composite_minimap2_index [file]   Path to composite minimap2 .mmi file to speed up analysis
            --keep_ref_id [str]                 String name of reference contig ID to keep during host removal (Default: 'MN908947.3')
            --min_length [i]                    Minimum length of fastq reads to keep (Default: 350)
            --max_length [i]                    Maximum length of fastq reads to keep (Default: 2400)
            --min_read_count [i]                Minimum read count required to output results (Default: 1)
            --flat                              Flag to flatten fastq output from named dirs to flat files in the output fastq_pass dir

        Example Command:
            nextflow run phac-nml/ncov-dehoster --human_ref ./hg38.fa -profile conda --nanopore --minimap2
                --fastq_directory ./fastqs 
                --run_name 'example_run_name'

                                        ------------------------------------------

    Nanostripper Pipeline Arguments:
      **NOTE: I DO NOT RECOMMEND RUNNING THIS PIPELINE AT THE MOMENT (NOT MAINTAINED), USE THE MINIMAP2 ONE ABOVE**
        Mandatory:
            --run_name [str]                Seperate results based on a run name input
            --fast5_directory [path]        Path to barcoded fast5 directories
            --nanostripper_tool_path [path] Path to nanostripper tool itself        
            --nanostripper_env_path [path]  Path to nanostripper conda environment (Has to be manually installed)
            --guppyCPU [path]               Path to guppy CPU conda environment to demultiplex (Proprietary, has to be manually installed)
            --guppyGPU [path]               Path to guppy GPU conda environment to basecall

        Optional:
            --min_length [i]                Minimum length of fastq reads to keep (Default: 350)
            --max_length [i]                Maximum length of fastq reads to keep (Default: 2400)
            --human_minimap2_index [file]   Path to minimap2 human indexed .mmi file to speed up analysis

        Example Command:
            nextflow run phac-nml/ncov-dehoster --human_ref ./hg38.fa -profile conda --nanopore --nanostripper
                --fast5_directory ./fast5_pass/
                --nanostripper_tool_path ./nanostripper
                --nanostripper_env_path ./conda/nanostripper_env
                --guppyCPU ./conda/guppyCPU_env/
                --guppyGPU ./conda/guppyGPU_env/
      
  """.stripIndent()
}
