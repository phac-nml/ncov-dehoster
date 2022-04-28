#!/bin/bash
set -eo pipefail

# Create a Cache Dir
mkdir -p conda_cache_dir
mkdir -p illumina_fastq

# Get Datasets
# Get partial human ref genome
wget https://raw.githubusercontent.com/DarianHole/test-datasets/master/partial_human_ref/partial_hg38_ref.fa.gz
gunzip partial_hg38_ref.fa.gz
# Small fastq paried files
wget https://raw.githubusercontent.com/DarianHole/test-datasets/master/illumina_fastq/illumina-18_S15_R1.fastq.gz https://raw.githubusercontent.com/DarianHole/test-datasets/master/illumina_fastq/illumina-18_S15_R2.fastq.gz -P ./illumina_fastq
wget https://raw.githubusercontent.com/DarianHole/test-datasets/master/illumina_fastq/illumina-4_S12_R1.fastq.gz https://raw.githubusercontent.com/DarianHole/test-datasets/master/illumina_fastq/illumina-4_S12_R2.fastq.gz -P ./illumina_fastq

# Run Flat Pipeline
rm -f /usr/share/miniconda/pkgs/cache/*.json # workaround for mamba-org/mamba#488

nextflow run ./main.nf \
    -profile conda,test \
    --cache ./conda_cache_dir \
    --illumina \
    --directory $PWD/illumina_fastq \
    --human_ref $PWD/partial_hg38_ref.fa

# Reset and Track
mv .nextflow.log artifacts/illumina.nextflow.log
rm -rf results work/ .nextflow*

echo "Done"
