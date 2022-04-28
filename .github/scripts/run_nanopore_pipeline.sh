#!/bin/bash
set -eo pipefail

# Create a Cache Dir
mkdir -p conda_cache_dir
mkdir -p nanopore_fastq

# Get Datasets
# Get partial human ref genome
wget https://raw.githubusercontent.com/DarianHole/test-datasets/master/partial_human_ref/partial_hg38_ref.fa.gz
gunzip partial_hg38_ref.fa.gz
# Small fastq files
wget https://raw.githubusercontent.com/DarianHole/test-datasets/master/nanopore_fastq/nanopore-1.fastq.gz -P ./nanopore_fastq
wget https://raw.githubusercontent.com/DarianHole/test-datasets/master/nanopore_fastq/nanopore-2.fastq.gz -P ./nanopore_fastq

# Run Flat Pipeline
rm -f /usr/share/miniconda/pkgs/cache/*.json # workaround for mamba-org/mamba#488

nextflow run ./main.nf \
    -profile conda,test \
    --cache ./conda_cache_dir \
    --nanopore \
    --minimap2 \
    --fastq_directory $PWD/nanopore_fastq \
    --run_name 'test-1-minimap2-flat' \
    --flat \
    --human_ref $PWD/partial_hg38_ref.fa

# Reset and Track
mv .nextflow.log artifacts/minimap2_flat.nextflow.log
rm -rf results work/ .nextflow*

# Run non-flat and test for cache dir working
nextflow run ./main.nf \
    -profile conda,test \
    --cache ./conda_cache_dir \
    --nanopore \
    --minimap2 \
    --fastq_directory $PWD/nanopore_fastq \
    --run_name 'test-2-minimap2' \
    --human_ref $PWD/partial_hg38_ref.fa

# Reset and Track
mv .nextflow.log artifacts/minimap2_expanded.nextflow.log
rm -rf results work/ .nextflow*

echo "Done"
