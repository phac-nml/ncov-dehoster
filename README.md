# ncov-dehoster

## About:
Nextflow pipeline that removes human reads from SARS-CoV-2 Illumina or Nanopore sequencing data. Basic details for the three available pipelines are as follows:

**Illumina Paired Pipeline \(v0.5.1)** - Competitive mapping approach using [bwa mem](http://bio-bwa.sourceforge.net/bwa.shtml) to remove human reads from input fastq file pairs while maintaining as many viral reads as possible.

**Nanopore Minimap2 Fastq Pipeline \(v0.5.1)** - Competitive mapping approach using [minimap2](https://github.com/lh3/minimap2) to remove human reads from either input fastq files or barcoded fastq directories while maintaining as many viral reads as possible. 
- Strict demultiplexing is *highly recommended* before running (unable to do so after and it improves downstream analyses)
- Optional fast5 dehosting available with argument `--fast5_directory [dir]`

**Nanopore Nanostripper Fast5 Pipeline \(v0.1.0-Depreciated)** - Dual mapping approach based around [nanostripper](https://github.com/nodrogluap/nanostripper) and [guppy](https://nanoporetech.com/nanopore-sequencing-data-analysis) to generate dehosted, demultiplexed fast5 and fastq files from input fast5 files. 

- *Currently Developmental*, the nanopore nanostripper dehosting pipeline is slow (full 96 sample run takes 6+ hours) and not setup well for external users which makes it **unrecommended** to try. You will have to provide a path to the guppy environment(s), the nanostripper environment and the nanostripper tool itself with the following arguments to generate a completely dehosted run and then it still may not work correctly:
    - `--guppyGPU <path/to/guppyGPU/env>`
    - `--guppyCPU <path/to/guppyCPU/env>`
    - `--nanostripper_env_path <path/to/nanostripper/env`
    - `--nanostripper_tool_path <path/to/nanostripper/tool`

------

## Index

[About](#about)

[Changelog](#changelog-highlights)

[Dependencies and Data](#dependencies-and-data-required)
- [Dependencies](#dependencies)
- [Data Required](#data-required)

[Quick Start](#quick-start)
- [Illumina](#illumina)
- [Nanopore Minimap2 Fastq](#minimap2-fastq-pipeline)
- [Nanopore Nanostripper Fast5 (Developmental)](#nanostripper-fast5-pipeline-developmental)

[Pipeline Full Details](#pipeline-details)
- [Illumina Paired Fastq Dehosting Pipeline](#illumina-paired-fastq-dehosting-pipeline)
- [Nanopore Fastq Minimap2 Dehosting and Regeneration Pipeline](#nanopore-fastq-minimap2-dehosting-and-regeneration-pipeline)
- [Nanopore Fast5 Nanostripper Dehosting and Regeneration Pipeline (Developmental)](#nanopore-fast5-nanostripper-dehosting-and-regeneration-pipeline-developmental)

[Downsampling](#downsampling)

[Profiles](#profiles)

[Other Notes](#other-notes)

------

## Changelog Highlights

#### Release v0.5.1
- Final Illumina fastq files now sorted by read names (QNAME field) rather than the chromosomal coordinates (default)
    - Addresses rare instance of output reads not being properly paired

#### Release v0.5.0
- Adjusted resource labels to match nf-core's along with adding in retries for specific processes

#### Release v0.4.0
- Added in optional downsampling with `--downsample` and a set of optional parameters
    - Random downsampling of final fastq reads with seqtk
    - Amplicon based downsampling of bam file with samtools
- Mamba added as its own profile instead of forcing it with conda

#### Release v0.3.0
- Adjusted how conda is implemented to support newer nextflow versions
- Fixed integration tests

#### Release v0.2.0
- Better user parameter options
    - Added `--keep_ref_id <ID>` to allow user to pick what reference ID they want to keep in output fastq files
- Process and tool version output added

#### Release v0.1.0
- Initial release

------

## Dependencies and Data Required

### Dependencies

#### Illumina Pipeline
- csvtk
- bwa
- pysam
- samtools
- seqtk

#### Nanopore Minimap2 Fastq Pipeline
- artic
- csvtk
- minimap2
- ont-fast5-api
- pip
- pysam
- samtools
- seqkt

#### Nanopore Nanostripper Fast5 Pipeline
- artic
- csvtk
- guppy-cpu
- guppy-gpu
- nanostripper

### Data Required

- Human Reference Sequence Fasta file
    - Get the reference genome from NCBI with the following:
    ```
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
    gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
    ```

    - Note that the newly released human reference genome (from April 2022 I believe) has not yet been tested but should work with the process

------

## Quick Start:

Simple commands to start dehosting along with the dependencies required . Greater detail on what each pipeline entails are found below in the individual [Pipeline Details Sections](#pipeline-details). You can also do `nextflow run phac-nml/ncov-dehoster --help` to get more info on the available processes

### Illumina

Run the *Illumina* pipeline with the following command:

```
nextflow run phac-nml/ncov-dehoster -profile conda --illumina --directory <path/to/paired_reads/dir> --human_ref <path/to/reference>
```

Minimum computational specs required for the Illumina dehosting pipeline are 3+ cpu and 16+G memory. More resources may be needed depending on the size of the input fastq files.

### Nanopore

#### Minimap2 Fastq Pipeline

Run the basic *Minimap2 Fastq* pipeline with the following command:

```
nextflow run phac-nml/ncov-dehoster -profile conda --nanopore --minimap2 --fastq_directory <path/to/fastqs> --human_ref <path/to/reference> --run_name 'whatever_you_want'
```

Minimum computational specs required for the minimap2 nanopore fastq host removal pipeline are 2+ cpu and 12+G memory. More resources may be needed depending on the size of the input fastq files.

#### Nanostripper Fast5 Pipeline (Depreciated)

*Requires conda* due to use of how the guppy environments are input as parameters at the moment (and that they cannot be made with conda as they are proprietary)

Run the full *Nanopore Nanostripper* pipeline with the following command:

```
nextflow run phac-nml/ncov-dehoster -profile conda --nanopore --nanostripper --fast5_directory <path/to/fast5_pass/> --human_ref <path/to/reference> --run_name 'whatever_you_want' --guppyCPU </path/to/conda_env/guppy-cpu/> --guppyGPU </path/to/conda_env/guppy-gpu/>
```

You can also just generate dehosted fast5 files with nanostripper with no guppy environment specified. Guppy is proprietary software of ONT Technologies so you must create your own environment for it

Minimum computational specs required for the Nanopore dehosting pipeline are 16+ cpu and 68+G memory. Configuring a GPU is highly recommended to run basecalling or the pipeline will take >12 hours to run the basecalling step (with a GPU it will be about 1-4 hours).

------

## Pipeline Details:

Full details for each pipeline including:
- Inputs
- Outputs
- Running
- Process
- Extra info

### Illumina Paired Fastq Dehosting Pipeline 

#### **Inputs**

Specify running Illumina dehosting with the `--illumina` flag.

Other arguments include:

| Parameter | Description | Default | Optional |
|-|-|-|-|
| directory | Directory containing paired fastq reads | None | No |
| human_ref | Fasta formatted human reference genome | None | No |
| cov2019_ref | Fasta formatted Sars-CoV-2 reference genome | data/nCoV-2019.reference.fasta | No |
| keep_ref_id | Reference genome header id to keep | MN908947.3 | No |
| keep_min_map_quality | Minimum mapping quality of covid reads to keep | 60 | No |
| remove_min_map_quality | Minimum mapping quality of the human reads to remove | 0 | No |
| composite_bwa_index | Directory containing BWA indexes for a composite human/viral reference --  **Speeds up analysis if given | None | Yes |
| downsample | Turn on downsampling of reads with seqtk | False | Yes |
| downsample_count | Approximate number of reads to be kept for each sample | 200,000 | Yes |
| downsample_seed | Integer seed to use for downsampling | 101 | Yes |
| downsample_amplicons | Utilize bed file to downsample sequencing amplicons with samtools | None | Yes |
| max_memory | Maximum memory to allow for a process | 256.GB | Yes |
| max_cpus | Maximum number of CPUs to allow for a process | 16 | Yes |
| max_time | Maximum time to allow a process to run for | 120.h | Yes |

#### **Running**

Full instructions on how to easily install and run the Illumina dehosting pipeline.

1. Setup all necessary resources:

    - [Conda](https://conda.io/en/latest/miniconda.html) with nextflow installed into an environment

    - A copy of the hg38 human reference genome

    - Folder with paired `fastq` or `fastq.gz` files to dehost

2. Activate the conda environment and run the pipeline

    ```
    nextflow run phac-nml/ncov-dehoster -profile conda --illumina --directory <path/to/reads> --human_ref <path/to/hg38.fa>
    ```

3. All subsequent runs can be extremely sped up using the `full/path/to/results/humanBWAIndex` folder as follows

    ```
    nextflow run phac-nml/ncov-dehoster -profile conda --illumina --directory <path/to/reads> --human_ref <path/to/hg38.fa> --composite_bwa_index </full/path/to/results/humanBWAIndex/>
    ```

**Note**: At the moment, you need to specify the full path to the made index or the pipeline will error out! But adding this removes the > 1 hour indexing step!

#### **Outputs**

Found in `./results/` directory, the outputs for the Illumina pipeline include:

- compositeMAPs --> Folder containing the composite mapping sorted BAM files along with each samples flagstats output

- dehostedBAMs --> Folder containing the dehosted composite BAM files

- dehosted_paired_fastqs --> Folder containing the final dehosted, paired fastq reads (main and final output)

- removal_summary --> CSV file containing read removal metrics

- process_versions --> YML file containing process names with their tools and versions associated for tracking

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
        - Sorted BAM file
        - flagstats output

    - Tools:
        - bwa mem
        - samtools sort

3. Host read removal and quality filtering

    Use the viral reference contig name (default MN908947.3) to determine which reads mapped to the human reference with >= the mapping quality score (default 0). Then pull out only proper pairs matching the viral reference and make sure they are above the quality wanted (default 60).

    - Inputs:
        - Sorted BAM file

    - Outputs:
        - Dehosted BAM file

    - Tools:
        - dehost_illumina.py
            - pysam

4. Regenerate fastq files

    Regenerate paired fastq files with samtools fastq

    - Inputs:
        - Dehosted BAM file

    - Outputs:
        - Dehosted, paired fastq reads files
    
    - Tools:
        - samtools fastq

#### **Additional Info**

- Using `--composite_bwa_index </path/to/composite_index_folder/` will speed up the analysis as the index won't need to be generated
    - Example generation: `cat human.fa sarsCoV2.fa > composite.fa && bwa index -a bwtsw composite.fa`

- If you are using a different SarsCoV2 reference genome than the one provided, make sure to match the `--keep_ref_id` so that the competitive mapping will pull out viral reads

- The `keep_min_map_quality` and `remove_min_map_quality` can be adjusted depending on the tolerance for human reads. The defaults are set at the most strict that then can be

- Host removal status can be checked with Kraken2 using a command similar to the following one each sample in the run:
    ```
    kraken2 --confidence 0.1 --db PATH/TO/kraken2_covid19_human_db/ --threads <THREADS> --report <SAMPLE>-REPORT.tsv --output <SAMPLE>-kraken.tsv <FASTQ_FILE_IN>
    ```

    - And then checking the total human reads with:
        ```
        grep -P "\s+9606\s+" -H *-kraken.tsv > human.tsv | wc -l human.tsv
        ```

------

### Nanopore Fastq Minimap2 Dehosting and Regeneration Pipeline

#### **Inputs**

Specify running Nanopore fastq dehosting using the minimap2 pipeline with both the --nanopore and --minimap2 flags as input.

Other arguments include:

| Parameter | Description | Default | Optional |
|-|-|-|-|
| fastq_directory | Directory containing either individual fastq files or barcoded fastq files | None | No |
| human_ref | Fasta formatted human reference genome | None | Yes - Need to pass --human_minimap2_index instead |
| composite_minimap2_index | Path to human minimap2 index (*.mmi) | None | Yes - Need to pass --human_ref instead |
| cov2019_ref | Fasta formatted Sars-CoV-2 reference genome | data/nCoV-2019.reference.fasta | No |
| run_name | Specify output run name to separate results | None | No |
| min_length | Minimum fastq read length to keep | 400 | Yes |
| max_length | Maximum fastq read length to keep | 2400 | Yes |
| min_read_count | Minimum read count required to output results | 1 | Yes
| fast5_directory | Directory of run associated fast5 files to be dehosted | None | Yes |
| keep_ref_id | Reference genome header id to keep | MN908947.3 | No |
| flat | Output flat fastq_pass folder instead of sample name subdirectories (better for folder input of named files) | False | Yes
| downsample | Turn on downsampling of reads with seqtk | False | Yes |
| downsample_count | Approximate number of reads to be kept for each sample | 100,000 | Yes |
| downsample_seed | Integer seed to use for downsampling | 101 | Yes |
| downsample_amplicons | Utilize bed file to downsample sequencing amplicons with samtools | None | Yes |

#### **Running**

Full instructions on how to easily install and run the Nanopore Minimap2 fastq dehosting pipeline.

1. Setup all necessary resources:

    - [Conda](https://conda.io/en/latest/miniconda.html) with nextflow installed into an environment OR all tools installed and available to skip using conda

    - A copy of the hg38 human reference genome OR a composite minimap2 index (mmi file) for the human reference genome and Sars-CoV-2 genome

    - One of the following as an input to dehost:
        - A directory with barcoded directories containing fastq files
        - A directory with named singular fastq files (one fastq file per sample)

2. (OPTIONAL) Strictly demultiplex the data (if using barcoded directories as input) before running the pipeline

    - Demultiplex the data using [guppy_barcoder](https://community.nanoporetech.com) (proprietary from ONT) and the `--require_barcodes_both_ends` flag to make sure that barcodes are present on each end of the reads.
        - There are some spots in the genome where this helps in resolution when analyzing the data
    
    - **Note** that the output fastq files are not able to be demultiplexed anymore so if you wish to strictly demultiplex the data it is advised to do it before running this host-removal pipeline.
        - May be able to be fixed in a later version

3. Ensure nextflow is available and then run the pipeline

    Basic command with conda to generate just dehosted fastq files out:
    ```
    nextflow run phac-nml/ncov-dehoster -profile conda --nanopore --minimap2 --fastq_directory <path/to/reads> --human_ref <path/to/hg38.fa>
    ```

    Command to generate dehosted fastq and fast5 files along with create a new sequencing_summary file:
    ```
    nextflow run phac-nml/ncov-dehoster -profile conda --nanopore --minimap2 --fastq_directory <path/to/reads> --human_ref <path/to/hg38.fa> --fast5_pass <path/to/fast5>
    ```

4. All subsequent runs can be sped up using the `full/path/to/results/compositeMinimapIndex/composite_ref.mmi` file as follows:

    ```
    nextflow run phac-nml/ncov-dehoster -profile conda --nanopore --minimap2 --fastq_directory <path/to/reads> --human_ref <path/to/hg38.fa> --composite_minimap2_index <full/path/to/results/compositeMinimapIndex/composite_ref.mmi>
    ```

**Note**: At the moment, you need to specify the full path to the made index or the pipeline will error out! But adding this removes the long indexing step!

#### **Outputs**

Found in `./results/<run_name>/run/` directory, the outputs for the Nanopore pipeline include:

- fastq_pass --> Folder containing the finished dehosted fastq files that pass the minimum read check (default: 1) which can be used for any other analyses
    - Separated by either barcode## or sample name depending upon input fastq folder structure
    - If `--flat` is passed, all fastq files will be in the main `fastq_pass` directory instead of named subdirectories.

- fast5_pass --> Folder containing the finished dehosted fast5 files that can be used for any other analyses
    - Separated by either barcode## or sample name depending upon input fastq folder structure
    - Only if `--fast5_directory` arg is passed

- sequencing_summary --> Simplified dehosted sequencing summary output only containing the read name and its specific fast5 file to allow data to be re-ran.
    - Only if `--fast5_directory` arg is passed

- removal_summary --> CSV file containing read removal metrics. Found in `./results/<run_name>/` instead
    - Shows the number and percentage of reads kept
    - Shows if the sample failed the minimum read filter

- process_versions --> YML file containing process names with their tools and versions associated for tracking

The output structure is setup as such so that the `run_name` organizes the sequencing data in the `run` folder and all analyses can be done in the `run_name` folder to separate out runs better.

#### **Process**

1. Composite genome generation and indexing with Minimap2

    Generate a composite reference genome and index it with minimap2 (skipped if --composite_minimap2_index is given)

    - Inputs:
        - Sars-CoV-2 reference fasta genome
        - Human reference fasta genome

    - Outputs:
        - composite reference index
    
    - Tools:
        - minimap2

2. Fastq size selection

    Filtering out fastq reads that do not fit the expected size range given (default 400 - 2400)

    - Inputs:
        - Fastq files or barcoded directories (based on what is given on the CL)

    - Outputs:
        - Size selected fastq file associated with its sample name

    - Tools:
        - artic guppyplex

3. Competitive Mapping

    Use competitive mapping to find human and chimeric reads while retaining known viral reads that may also map weakly to the human reference.

    - Inputs:
        - Size selected fastq file and its sample name
        - Composite minimap2 reference (*.mmi file)

    - Outputs:
        - Sorted BAM file

    - Tools:
        - minimap2
        - samtools

4. Host read removal and quality filtering

    Use the viral reference contig name (default MN908947.3) to determine which reads mapped to the human reference with >= the mapping quality score (default 0). Then pull out only reads mapping to the viral reference and make sure they are above the quality wanted (default 60).

    - Inputs:
        - Sorted BAM file

    - Outputs:
        - Dehosted BAM file
        - Sample specific removal CSV

    - Tools:
        - dehost_nanopore.py
            - pysam

5. Regenerate fastq files

    Regenerate paired fastq files with samtools fastq

    - Inputs:
        - Dehosted BAM file

    - Outputs:
        - Dehosted fastq file
    
    - Tools:
        - samtools fastq

6. (OPTIONAL) Regenerate fast5 files and the sequencing summary file if `--fast5_directory` argument is given

    Regenerate dehosted fast5 files from the fastq files using fast5_subset

    - Inputs:
        - Dehosted fastq files
        - Fast5 directory associated with the run

    - Outputs:
        - dehosted fast5 files
        - sequencing_summary.txt
    
    - Tools
        - awk
        - fast5_subset

#### **Additional Info**

- Providing the composite minimap2 index genereated in a previous analysis will speed other ones

- Inputs can be fastq barcode directories or a directory of fastq files
    - Example barcode directory:
        ```
        fastq_in
        ├── barcode01
        │   ├── X_pass_barcode01_556763c2_0.fastq
        │   ├── X_pass_barcode01_556763c2_10.fastq
        │   └── X_pass_barcode01_556763c2_11.fastq
        ├── barcode02
        │   ├── X_pass_barcode02_556763c2_0.fastq
        │   ├── X_pass_barcode02_556763c2_10.fastq
        │   └── X_pass_barcode02_556763c2_11.fastq
        etc.
        ```
    - Example fastq directory:
        ```
        fastq_in
        ├── ABC.fastq
        ├── DEF.fastq
        ├── GHI.fastq
        ├── JKL.fastq
        etc.
        ```
- Host removal status can be checked with Kraken2 using a command similar to the following one each sample in the run:
    ```
    kraken2 --confidence 0.1 --db PATH/TO/kraken2_covid19_human_db/ --threads <THREADS> --report <SAMPLE>-REPORT.tsv --output <SAMPLE>-kraken.tsv <FASTQ_FILE_IN>
    ```

    - And then checking the total human reads with:
        ```
        grep -P "\s+9606\s+" -H *-kraken.tsv > human.tsv | wc -l human.tsv
        ```

------

### Nanopore Fast5 Nanostripper Dehosting and Regeneration Pipeline (Developmental)

#### **Inputs**

*At the moment, the [minimap2 pipeline](#nanopore-fastq-minimap2-dehosting-and-regeneration-pipeline) is recommended due to its ease of use.* This one still works, it is just a fair bit slower and a lot harder to set up due to some of the inputs relying on the user to set them.

Specify running Nanopore dehosting with both the --nanopore and --nanostripper parameters as input.

Other parameters include:

| Parameter | Description | Default | Optional |
|-|-|-|-|
| fast5_directory | Directory containing barcoded fast5 files (non-barcoded not yet supported) | None | No |
| human_ref | Fasta formatted human reference genome | None | Yes - Need to pass human_minimap2_index instead  |
| human_minimap2_index | Path to human minimap2 index (*.mmi) | None | Yes - Need human_ref |
| cov2019_ref | Fasta formatted Sars-CoV-2 reference genome | data/nCoV-2019.reference.fasta | No |
| run_name | Specify output run name to separate results | None | No |
| min_length | Minimum fastq read length (based on scheme) | 400 | Yes |
| max_length | Maximum fastq read length (based on scheme) | 2400 | Yes |
| guppy_cpu | Path to Guppy CPU conda environment | None | Yes - Will only de-host fast5 files |
| guppy_gpu | Path to Guppy GPU environment (to speed up basecalling only) | None | Yes - Also need to give CPU env |

#### **Outputs**

Found in `./results/<run_name>/run/` directory, the outputs for the Nanopore pipeline include:

- dehosted_fast5 --> Folder containing the dehosted fast5 results of nanostripper for each barcode

- fast5_pass --> Folder containing the finished dehosted demultiplexed fast5 files that can be used to run another analysis

- fastq_pass --> Folder containing the finished dehosted demultiplexed fastq files that can be used for another analysis

- sequencing_summary --> Simplified dehosted sequencing summary output only containing the read name and its specific fast5 file

- removal_summary --> CSV file containing read removal metrics. Found in `./results/<run_name>/` instead

The output structure is setup as such so that the `run_name` organizes the sequencing data in the `run` folder and all analyses can be done in the `run_name` folder to separate out runs better.

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

    Using guppy to re-basecall the remaining reads. This is done so that the most up to date version of Guppy is used and the best results are obtained as not all of the instraments are updated. A guppy environment must be provided as a parameter (at least CPU) to run this step and if none is given, the pipeline will not continue and you will only have the dehosted fast5 files.

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

## Downsampling

Added in v0.4.0, allows a user to either randomly downsample their final fastq files or downsample their dehosted bam files based on a given amplicon bed file.

Parameters are available in both the illumina and nanopore data pipelines and are as follows:
| Parameter | Description | Default | Optional |
|-|-|-|-|
| downsample | Turn on downsampling of reads with seqtk | False | Yes |
| downsample_count | Approximate number of reads to be kept for each sample | 200,000 Illumina - 100,000 Nanopore | Yes |
| downsample_seed | Integer seed to use for downsampling | 101 | Yes |
| downsample_amplicons | Utilize bed file to downsample sequencing amplicons with samtools | None | Yes |

To use downsampling, turn it on by passing in `--downsample` to your command and then adjust the number of reads and the seed with the other available args. Downsampling is done using [seqtk](https://github.com/lh3/seqtk) which randomly downsamples the fastq files based on a given seed.

Passing in the `--downsample_amplicons <BEDFILE>` will utilize a different approach to downsample the data. Instead of running on the dehosted fastq files, it is run on the dehosted bam file with samtools. Each definied region will be downsampled to a maximum read count based on the input `downsample_count` and the number of amplicons in the file which is calculated at run time. For example if `downsample_count=100000` and `amplicon_count=30`, then each amplicon would be downsampled to ~3334 reads.

The amplicon bed file should be structured as such:
```
MN908947.3	54	1205	1	nCoV-2019_1	+
MN908947.3	1128	2266	2	nCoV-2019_2	+
MN908947.3	2179	3257	3	nCoV-2019_1	+
MN908947.3	3166	4262	4	nCoV-2019_2	+
ContigID	AmpStart	AmpEnd	Amp#	Pool	Direction
```

Although only the first 3 columns are used to downsample so the remaining ones can be excluded if wanted.

------

## Profiles

Profiles are a set of configuration attributes that can be activated (as many as you want as long as they are comma separated) when launching a pipeline. More information can be found [here](https://www.nextflow.io/docs/latest/config.html?highlight=profile#config-profiles)

If no profile is given, the Illumina and nanopore minimap2 pipelines will still function correctly and will run all processes locally so long as the system has all of the dependencies installed and meets the computational requirements. The nanopore nanostripper pipeline currently requires a guppy conda environment to re-basecall and demultiplex the fast5 files.

Profiles can be activated with `-profile <profile_name>`. Below are the profiles currently available with this pipeline

### Conda

The conda profile can be activated by passing `-profile conda`.

This profile uses conda environments that are created on the fly or found in the given cache directory through the cache parameter (`--cache <path/to/conda/envs/>`) to control the dependencies and run the analysis.

If the environments are having trouble being made, it is recommended to use mamba to install them

### Mamba

The conda profile can be activated by passing `-profile mamba`.

This profile uses conda environments that are created on the fly with mamba or found in the given cache directory through the cache parameter (`--cache <path/to/conda/envs/>`) to control the dependencies and run the analysis.

### Custom

Custom profiles are specified in the main [nextflow config file](https://github.com/phac-nml/ncov-dehoster/blob/d553cccb64df0218deadc0f455c819f5469498c8/nextflow.config#L49) and allow for more specific allocation of resources along with other executors to be used (ex. slurm)

ex. [NML](https://github.com/phac-nml/ncov-dehoster/blob/master/conf/custom/nml.config) config uses slurm with specific options

**Additional configs can be added fairly easily by request or on your own fork of the repo

------

## Other Notes:

A fair amount of testing has gone into the pipeline but there may still be bugs or other issues. If you discover any please make an issue.
