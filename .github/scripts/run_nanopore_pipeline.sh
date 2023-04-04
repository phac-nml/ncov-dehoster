#!/bin/bash
set -eo pipefail

# Create a Cache Dir
mkdir -p conda_cache_dir

### Run Flat Pipeline ###
nextflow run ./main.nf \
    -profile conda,test \
    --cache ./conda_cache_dir \
    --nanopore \
    --minimap2 \
    --fastq_directory $PWD/.github/data/nanopore/ \
    --run_name 'test-1-minimap2-flat' \
    --flat \
    --composite_minimap2_index $PWD/.github/data/minimap2_index/sars-cov-2_ref.mmi

### Check Outputs ###
# 1. Num Poor Reads
READS=`awk -F, '$1 == "nanopore-1" {print $3}' ./results/test-1-minimap2-flat/removal_summary.csv`
if [[ "$READS" != "2" ]]; then 
    echo "Incorrect output: Number of poor reads found"
    echo "  Expected: 2, Got: $READS"
    exit 1
fi

# 2. Total Reads Kept
READS=`awk -F, '$1 == "nanopore-2" {print $4}' ./results/test-1-minimap2-flat/removal_summary.csv`
if [[ "$READS" != "2013" ]]; then 
    echo "Incorrect output: Number of human reads found"
    echo "  Expected: 2013, Got: $READS"
    exit 1
fi

# Reset and Track
mv .nextflow.log artifacts/minimap2_flat.nextflow.log
rm -rf results work/ .nextflow*

### Run non-flat and test for cache dir working and human ref working ###
nextflow run ./main.nf \
    -profile conda,test \
    --cache ./conda_cache_dir \
    --nanopore \
    --minimap2 \
    --fastq_directory $PWD/.github/data/nanopore/ \
    --run_name 'test-2-minimap2-partial-human' \
    --human_ref $PWD/.github/data/partial_hg38_ref.fa

### Check Outputs ###
# 1. Num Poor Reads
READS=`awk -F, '$1 == "nanopore-1" {print $3}' ./results/test-2-minimap2-partial-human/removal_summary.csv`
if [[ "$READS" != "2" ]]; then 
    echo "Incorrect output: Number of poor reads found"
    echo "  Expected: 2, Got: $READS"
    exit 1
fi

# 2. Total Reads Kept
READS=`awk -F, '$1 == "nanopore-2" {print $4}' ./results/test-2-minimap2-partial-human/removal_summary.csv`
if [[ "$READS" != "2013" ]]; then 
    echo "Incorrect output: Number of human reads found"
    echo "  Expected: 2013, Got: $READS"
    exit 1
fi

# Reset and Track
mv .nextflow.log artifacts/minimap2_expanded.nextflow.log
rm -rf results work/ .nextflow*

echo "Done"
