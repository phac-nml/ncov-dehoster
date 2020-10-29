# ncov-dehoster

## About:
Nextflow pipeline that removes human reads from SARS-CoV-2 Illumina and Nanopore sequencing data.

**Illumina** - Competitive mapping approach with bwa mem to remove human reads from the input fastq files while maintaining as many viral reads as possible

**Nanopore** - Dual mapping approach based around [nanostripper](https://github.com/nodrogluap/nanostripper) and then subsequently [guppy](https://nanoporetech.com/nanopore-sequencing-data-analysis) to generate de-hosted, demultiplexed fast5 and fastq files from input fast5 files

------

## Quick Start:

Simple commands to start de-hosting. Greater detail on what each entails is found below in the Pipeline Details sections

### Illumina

Run the **Illumina** pipeline with the following command:

```
nextflow run phac-nml/ncov-dehoster -profile conda --illumina --directory <path/to/paired_reads/dir> --human_ref <path/to/reference>
```

### Nanopore

*Requires conda* due to use of guppy environments as parameters

Run the full **Nanopore** pipeline with the following command:

```
nextflow run phac-nml/ncov-dehoster -profile conda --nanopore --directory <path/to/fast5_pass/> --human_ref <path/to/reference> --run_name 'whatever_you_want' --guppyCPU </path/to/conda_env/guppy-4.0.11-cpu/>
```

You can also just generate de-hosted fast5 files with nanostripper with no guppy environment specified. Guppy is proprietary software of ONT Technologies so you must create you're own environment for it

------

## Pipeline Details:

### Illumina Paired Fastq De-hosting Pipeline 

#### **Inputs**

Specify running Illumina dehosting with the `--illumina` parameter.

Other parameters include:

| Parameter | Description | Default | Optional |
|-|-|-|-|
| directory | Directory containing paired fastq reads | None | No |
| human_ref | Fasta formatted human reference genome | None | No |
| cov2019_ref | Fasta formatted Sars-CoV-2 reference genome | data/nCoV-2019.reference.fasta | No |
| covid_ref_id | Sars-CoV-2 reference genome header id | MN908947.3 | No |
| keep_min_map_quality | Minimum mapping quality of covid reads to keep | 60 | No |
| remove_min_map_quality | Minimum mapping quality of the human reads to remove | 0 | No |
| composite_bwa_index | Directory containing BWA indexes for a composite human/viral reference --  **Speeds up analysis if given | None | Yes |

#### **Outputs**

Found in `./results/` directory, the outputs for the Illumina pipeline include:

- compositeMAPs --> Folder containing the composite mapping sorted bam files along with each samples flagstats output

- dehostedBAMs --> Folder containing the de-hosted composite bam files

- dehosted_paired_fastqs --> Folder containing the final de-hosted, paired fastq reads (main and final output)

- removal_summary --> CSV file containing read removal metrics (with more to come)

#### **Process**

1. Composite genome generation and Indexing

    Generate bwa index to allow competitive mapping with bwa mem.

    - Inputs:
        - Sars-CoV-2 reference fasta genome
        - Human reference fasta genome

    - Outputs:
        - Folder of bwa mem indexes of the composite genome
    
    - Tools:
        - bwa index

2. Competitive mapping

    Use competitive mapping to find human and chimeric reads while retaining known viral reads that may also map weakly to the human reference.

    - Inputs:
        - Bwa mem indexes (`.amb`, `.ann`, `.pac`, `.sa`)
        - Paired fastq reads

    - Outputs:
        - Sorted bam file
        - flagstats output

    - Tools:
        - bwa mem
        - samtools sort

3. Host read removal and quality filtering

    Use the viral reference contig name (default MN908947.3) to determine which reads mapped to the human reference with >= the mapping quality score (default 0). Then pull out only proper pairs matching the viral reference and make sure they are above the quality wanted (default 60).

    - Inputs:
        - Sorted bam file

    - Outputs:
        - Dehosted bam file

    - Tools:
        - dehost.py

4. Regenerate fastq files

    Regenerate paired fastq files with samtools fastq

    - Inputs:
        - Dehosted bam file

    - Outputs:
        - Dehosted, paired fastq reads files
    
    - Tools:
        - samtools fastq

#### **Additional Info**

- Using `--composite_bwa_index </path/to/composite_index_folder/` will speed up the analysis as the index won't need to be generated
    - Example generation: `cat human.fa sarsCoV2.fa > composite.fa && bwa index -a bwtsw composite.fa`

- If you are using a different SarsCoV2 reference genome than the one provided, make sure to match the `--covid_ref_id` so that the competitive mapping will pull out viral reads

- The `keep_min_map_quality` and `remove_min_map_quality` can be adjusted depending on the tolerance for human reads. The defaults are set at the most strict that then can be

------

### Nanopore De-hosting Pipeline

#### **Inputs**

| Parameter | Description | Default | Optional |
|-|-|-|-|
| directory | Directory containing barcoded fast5 files (non-barcoded not yet supported) | None | No |
| human_ref | Fasta formatted human reference genome | None | Yes - Need human_minimap2_index |
| human_minimap2_index | Path to human minimap2 index (*.mmi) | None | Yes - Need human_ref |
| cov2019_ref | Fasta formatted Sars-CoV-2 reference genome | data/nCoV-2019.reference.fasta | No |
| min_length | Minimum fastq read length (based on scheme) | 1600 | No |
| max_length | Maximum fastq read length (based on scheme) | 2400 | No |
| run_name | Specify output run name to separate results | None | No |
| guppy_cpu | Path to Guppy CPU conda environment | None | Yes - Will only de-host fast5 files |
| guppy_gpu | Path to Guppy GPU environment (to speed up basecalling only) | None | Yes - Also need to give CPU env |

#### **Outputs**

Found in `./results/<run_name>/` directory, the outputs for the Nanopore pipeline include:

- dehosted_fast5 --> Folder containing the dehosted fast5 results of nanostripper for each barcode

- fast5_pass --> Folder containing the finished de-hosted demultiplexed fast5 files that can be used to run another analysis

- fastq_pass --> Folder containing the finished de-hosted demultiplexed fastq files that can be used for another analysis

- sequencing_summary --> Simplified dehosted sequencing summary output only containing the read name and its specific fast5 file

#### **Process**

1. Nanostipper

    Use nanostripper to remove human reads from the input fast5 files with one nanostripper instance running per barcode

    - Inputs:
        - Sars-CoV-2 reference fasta genome
        - Human reference fasta genome or Human minimap2 indexed genome (speeds it up as only need to index once)

    - Outputs:
        - Folder barcoded and dehosted fast5 files
    
    - Tools:
        - nanostripper
        - minimap2 (inside nanostripper)

2. Guppy basecalling

    Using guppy to re-basecall the remaining reads. This is done so that the most up to date version of Guppy is used and the best results are obtained as not all of the instraments are updated. A guppy environment must be provided as a parameter (at least CPU) to run this step and if none is given, the pipeline will not continue and you will only have the de-hosted fast5 files.

    - Inputs:
        - Dehosted barcoded fast5 files

    - Outputs:
        - Dehosted only fastq files

    - Tools:
        - guppy_basecaller CPU or GPU (faster)

3. Size selection

    Filtering out dehosted fastq reads that do not fit the expected size of the primer scheme. All of the dehosted fastq reads are combined into one massive fastq file to run artic guppy_plex on and filter out reads of the wrong size

    - Inputs:
        - All dehosted fastq files combined into one big fastq file

    - Outputs:
        - Size selected dehosted fastq file (one big one still)

    - Tools:
        - artic guppy_plex

4. demultiplex fastq files
    
    The filtered and re-basecalled fastq files have to be re-demultiplexed incase they are a different barcode than originally thought or to get some of the unclassified reads to a barcode with the improvements from guppy_basecaller.

    - Inputs:
        - Size selected dehosted fastq file

    - Outputs:
        - fastq_pass - Barcoded folders containing their respective size selected dehosted fastq files 

    - Tools:
        - guppy_barcoder

5. Regenerate and demultiplex fast5 files and the simplified sequencing summary

    Using the fully finished fastq files in the `fastq_pass` directory, re-demultiplex the fast5 files to match and generate a new simplified sequencing summary file for the run that just contains the read name and fast5 file its in. This is accomplished with fast5_subset.

    - Inputs:
        - fastq_pass folder with barcoded size selected dehosted fastq files
        - Dehosted only fast5 folders from nanostripper

    - Outputs:
        - fast5_pass - Barcoded dehosted fast5 files matching the fastq ones
        - sequencing_summary.txt

    - Tools:
        - fast5_subset

#### **Additional Info**

- `-profile conda` is required due to specifiying `guppy` environment(s) as a parameter(s)
    - if not given, pipeline exits after nanostripper with dehosted fast5 files only

- The `--guppyCPU` parameter is required for demultiplexing so at minimum that environment must be made

- We re-basecall to get the benefits of the new `guppy` changes which improve accuracy in the output fastq files

- The `--guppyGPU` parameter can be given as a way to speed up basecalling but it is not required and if not given, basecalling will be done with the specified `--guppyCPU` environment

------

### Profiles

Profiles are a set of configuration attributes that can be activated (as many as you want as long as they are comma separated) when launching a pipeline. More information can be found [here](https://www.nextflow.io/docs/latest/config.html?highlight=profile#config-profiles)

If no profile is given, the Illumina pipeline will still function correctly and will run locally so long as the system has all of the dependencies installed and meets the computational requirements for the processes. The nanopore pipeline currently requires a guppy conda environment to re-basecall and demultiplex the fast5 files

Profiles can be activated with `-profile <profile_name>`. Below are the profiles currently available with this pipeline

#### **Conda**

The conda profile can be activated by passing `-profile conda`.

This profile uses conda environments that are created on the fly or given by the cache parameter as `--cache <path/to/conda/envs/>` to control the dependencies and run the analysis

#### **Custom**

Custom profiles are specified in the main [nextflow config file](https://github.com/phac-nml/ncov-dehoster/blob/d553cccb64df0218deadc0f455c819f5469498c8/nextflow.config#L49) and allow for more specific allocation of resources along with other executors to be used (ex. slurm)

ex. [NML](https://github.com/phac-nml/ncov-dehoster/blob/master/conf/custom/nml.config) config uses slurm with specific options

**Additional configs can be added fairly easily by request

------

## Other Notes:

- This is an initial version and there are potentially still bugs and changes to be made, although it seems to be working as intended 

## Upcoming (Potentially) Additions and Changes:

- More info in `removal_summary.csv` output for Illumina data

- Removal info for Nanopore data

- Single-end fastq de-hosting (maybe)

- Non-barcoded fast5 de-hosting (maybe)

- No re-basecalling fastq generation (nanostripper only with no basecalling or demultiplexing)

- Single directory fast5 de-hosting (maybe)
