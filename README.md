# ncov-dehoster

## About:
Nextflow pipeline that removes human reads from SARS-CoV-2 Illumina and Nanopore sequencing data.

**Illumina** - Competitive mapping approach with bwa mem to remove human reads from the input fastq files while maintaining as many viral reads as possible

**Nanopore** - Dual mapping approach based around [nanostripper](https://github.com/nodrogluap/nanostripper) and then subsequently [guppy](https://nanoporetech.com/nanopore-sequencing-data-analysis) to generate de-hosted, de-multiplexed fast5 and fastq files

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

### Profiles

Profiles are a set of configuration attributes that can be activated (as many as you want as long as they are comma separated) when launching a pipeline. More information can be found [here](https://www.nextflow.io/docs/latest/config.html?highlight=profile#config-profiles)

If no profile is given, the Illumina pipeline will still function correctly and will run locally so long as the system has all of the dependencies installed and meets the computational requirements for the processes. Nanopore currently requires a guppy conda environment to re-basecall

Profiles can be activated with `-profile <profile_name>`. Below are the profiles currently available with this pipeline

#### **Conda**

The conda profile can be activated by passing `-profile conda`.

This profile uses conda environments that are created on the fly or given by the cache parameter as `--cache <path/to/conda/envs/>` to control the dependencies and run the analysis

#### **Custom**

Custom profiles are specified in the main [nextflow config file](https://github.com/phac-nml/ncov-dehoster/blob/d553cccb64df0218deadc0f455c819f5469498c8/nextflow.config#L49) and allow for more specific allocation of resources along with other executors to be used (ex. slurm)

ex. [NML](https://github.com/phac-nml/ncov-dehoster/blob/master/conf/custom/nml.config) config uses slurm with specific options

**Additional configs can be added fairly easily by request

------

### Illumina Paired De-hosting Pipeline 

#### **Inputs**

Specify running Illumina dehosting with the `--illumina` parameter.

Other parameters include:

| Parameter | Description | Default | Optional |
|-|-|-|-|
| directory | Directory containing paired fastq reads | None | No |
| human_ref | Fasta formatted human reference genome | None | No |
| cov2019_ref | Fasta formatted Sars CoV2 reference genome | data/nCoV-2019.reference.fasta | No |
| covid_ref_id | Sars CoV2 reference genome header id | MN908947.3 | No |
| keep_min_map_quality | Minimum mapping quality of covid reads to keep | 60 | No |
| remove_min_map_quality | Minimum mapping quality of the human reads to remove | 0 | No |
| composite_bwa_index | Directory containing BWA indexes for a composite human/viral reference --  **Speeds up analysis if given | None | Yes |

#### **Outputs**

Found in `./results/` the outputs for the Illumina pipeline include:

- compositeMAPs --> Folder containing the composite mapping sorted bam files along with each samples flagstats output

- dehostedBAMs --> Folder containing the de-hosted composite bam files

- dehosted_paired_fastqs --> Folder containing the final de-hosted, paired fastq reads (main and final output)

- removal_summary --> CSV file containing read removal metrics (with more to come)

#### **Process**

put some pathway picture and explanations

#### **Additional Info**

- Using `--composite_bwa_index </path/to/composite_index_folder/` will speed up the analysis as the index won't need to be generated
    - Example generation: `cat human.fa sarsCoV2.fa > composite.fa && bwa index -a bwtsw composite.fa`

- If you are using a different SarsCoV2 reference genome than the one provided, make sure to match the `--covid_ref_id` so that the competitive mapping will pull out viral reads

- The `keep_min_map_quality` and `remove_min_map_quality` can be adjusted depending on the tolerance for human reads. The defaults are set at the most strict that then can be

------

### Nanopore De-hosting Pipeline

#### **Inputs**

table soon

#### **Outputs**

- files are outputted! 

#### **Process**

put some pathway picture and explanations

#### **Additional Info**

something here

------

## Other Notes:

- This is an initial version and there are potentially still bugs and changes to be made, although it seems to be working as intended currently

## Upcoming (Hopefully) Additions and Changes:

- More info in `removal_summary.csv` output for Illumina data

- Removal info for Nanopore data

- Single-end fastq de-hosting (maybe)

- Single directory fast5 de-hosting (maybe)
