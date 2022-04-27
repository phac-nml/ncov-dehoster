#!/bin/bash
set -eo pipefail

# Create a Cache Dir
mkdir -p conda_cache_dir
mkdir -p nanopore_fastq

# Get Datasets
# Covid minimap2 index (not using composite ref due to size, will work without just for the testing)
wget https://raw.githubusercontent.com/DarianHole/test-datasets/master/minimap2/sars-cov-2_ref.mmi
# Small fastq files
wget https://raw.githubusercontent.com/DarianHole/test-datasets/master/nanopore_fastq/nanopore-1.fastq.gz -P ./nanopore_fastq
wget https://raw.githubusercontent.com/DarianHole/test-datasets/master/nanopore_fastq/nanopore-2.fastq.gz -P ./nanopore_fastq

# Run Flat Pipeline
nextflow run ./main.nf \
    -profile conda,test \
    --cache ./conda_cache_dir \
    --nanopore \
    --minimap2 \
    --fastq_directory $PWD/nanopore_fastq \
    --run_name 'test-1-minimap2-flat' \
    --flat \
    --composite_minimap2_index $PWD/sars-cov-2_ref.mmi

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
    --composite_minimap2_index $PWD/sars-cov-2_ref.mmi

# Reset and Track
mv .nextflow.log artifacts/minimap2_expanded.nextflow.log
rm -rf results work/ .nextflow*

echo "Done"
