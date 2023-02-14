#!/bin/bash
set -eo pipefail

# Create a Cache Dir
mkdir -p conda_cache_dir

# Workaround for mamba-org/mamba#488
rm -f /usr/share/miniconda/pkgs/cache/*.json

# Run Flat Pipeline
nextflow run ./main.nf \
    -profile conda,test \
    --cache ./conda_cache_dir \
    --nanopore \
    --minimap2 \
    --fastq_directory $PWD/.github/data/nanopore/ \
    --run_name 'test-1-minimap2-flat' \
    --flat \
    --composite_minimap2_index $PWD/.github/data/minimap2_index/sars-cov-2_ref.mmi

# Reset and Track
mv .nextflow.log artifacts/minimap2_flat.nextflow.log
rm -rf results work/ .nextflow*

# Run non-flat and test for cache dir working and human ref working
nextflow run ./main.nf \
    -profile conda,test \
    --cache ./conda_cache_dir \
    --nanopore \
    --minimap2 \
    --fastq_directory $PWD/.github/data/nanopore/ \
    --run_name 'test-2-minimap2-partial-human' \
    --human_ref $PWD/.github/data/partial_hg38_ref.fa

# Reset and Track
mv .nextflow.log artifacts/minimap2_expanded.nextflow.log
rm -rf results work/ .nextflow* partial_hg38_ref.fa

echo "Done"
