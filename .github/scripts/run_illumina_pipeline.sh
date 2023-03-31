#!/bin/bash
set -eo pipefail

# Create a Cache Dir
mkdir -p conda_cache_dir

### Run Illumina Pipeline ###
nextflow run ./main.nf \
    -profile conda,test \
    --cache ./conda_cache_dir \
    --illumina \
    --directory $PWD/.github/data/illumina/ \
    --human_ref $PWD/.github/data/partial_hg38_ref.fa

### Check Pipeline Outputs ###
# 1. Num Human Reads
READS=`awk -F, '$1 == "illumina-18_S15" {print $2}' ./results/removal_summary.csv`
if [[ "$READS" != "2" ]]; then 
    echo "Incorrect output: Number of human reads found"
    echo "  Expected: 2, Got: $READS"
    exit 1
fi

# 2. Total Reads Kept
READS=`awk -F, '$1 == "illumina-4_S12" {print $4}' ./results/removal_summary.csv`
if [[ "$READS" != "3980" ]]; then 
    echo "Incorrect output: Number of human reads found"
    echo "  Expected: 3980, Got: $READS"
    exit 1
fi

# Reset and Track
mv .nextflow.log artifacts/illumina.nextflow.log
rm -rf results work/ .nextflow*

echo "Done"
