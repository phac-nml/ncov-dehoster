#!/bin/bash
set -eo pipefail

# Create a Cache Dir
mkdir -p conda_cache_dir

# Workaround for mamba-org/mamba#488
rm -f /usr/share/miniconda/pkgs/cache/*.json

# Run Illumina Pipeline
nextflow run ./main.nf \
    -profile conda,test \
    --cache ./conda_cache_dir \
    --illumina \
    --directory $PWD/.github/data/illumina/ \
    --human_ref $PWD/.github/data/partial_hg38_ref.fa

# Reset and Track
mv .nextflow.log artifacts/illumina.nextflow.log
rm -rf results work/ .nextflow* partial_hg38_ref.fa

echo "Done"
