#!/bin/bash
set -eo pipefail

# Create a Cache Dir
mkdir -p conda_cache_dir

# --------------------------------------------------------------------- #
### Run Illumina Pipeline ###
nextflow run ./main.nf \
    -profile mamba,test \
    --cache ./conda_cache_dir \
    --illumina \
    --directory $PWD/.github/data/illumina/ \
    --human_ref $PWD/.github/data/partial_hg38_ref.fa

### Check Pipeline Outputs ###
# 1. Num Human Reads
READS=`awk -F, '$1 == "illumina-18_S15" {print $2}' ./results/removal_summary.csv`
if [[ "$READS" != "2" ]]; then 
    echo "Incorrect output: Number of Human Reads"
    echo "  Expected: 2, Got: $READS"
    exit 1
fi

# 2. Total Reads Kept
READS=`awk -F, '$1 == "illumina-4_S12" {print $4}' ./results/removal_summary.csv`
if [[ "$READS" != "3980" ]]; then 
    echo "Incorrect output: Number of Paired Reads Kept"
    echo "  Expected: 3980, Got: $READS"
    exit 1
fi

# Reset and Track
mv .nextflow.log artifacts/illumina.nextflow.log
rm -rf results work/ .nextflow*

# --------------------------------------------------------------------- #
### Run Illumina Pipeline with Downsampling ###
nextflow run ./main.nf \
    -profile mamba,test \
    --cache ./conda_cache_dir \
    --illumina \
    --directory $PWD/.github/data/illumina/ \
    --human_ref $PWD/.github/data/partial_hg38_ref.fa \
    --downsample \
    --downsample_count 100 \
    --downsample_seed 42

### Check Pipeline Outputs ###
# 1. Num Human Reads
READS=`awk -F, '$1 == "illumina-18_S15" {print $2}' ./results/removal_summary.csv`
if [[ "$READS" != "2" ]]; then 
    echo "Incorrect output: Number of Human Reads"
    echo "  Expected: 2, Got: $READS"
    exit 1
fi

# 2. Total Reads Kept
READS=`awk -F, '$1 == "illumina-4_S12" {print $4}' ./results/removal_summary.csv`
if [[ "$READS" != "100" ]]; then 
    echo "Incorrect output: Number of Paired Reads Kept"
    echo "  Expected: 100, Got: $READS"
    exit 1
fi

# 3. Downsample Maximum
READS=`awk -F, '$1 == "illumina-4_S12" {print $6}' ./results/removal_summary.csv`
if [[ "$READS" != "100" ]]; then 
    echo "Incorrect output: Downsampling Maximum"
    echo "  Expected: 100, Got: $READS"
    exit 1
fi

# 4. Downsample Seed
SEED=`awk -F, '$1 == "illumina-4_S12" {print $7}' ./results/removal_summary.csv`
if [[ "$SEED" != "42" ]]; then 
    echo "Incorrect output: Downsampling Seed"
    echo "  Expected: 42, Got: $SEED"
    exit 1
fi

# Reset and Track
mv .nextflow.log artifacts/illumina_downsampled.nextflow.log
rm -rf results work/ .nextflow*

# --------------------------------------------------------------------- #
### Run Illumina Pipeline with Downsampling using Amplicons ###
nextflow run ./main.nf \
    -profile mamba,test \
    --cache ./conda_cache_dir \
    --illumina \
    --directory $PWD/.github/data/illumina/ \
    --human_ref $PWD/.github/data/partial_hg38_ref.fa \
    --downsample \
    --downsample_count 100 \
    --downsample_seed 42 \
    --downsample_amplicons $PWD/.github/data/amplicon.bed

### Check Pipeline Outputs ###
# 1. Num Human Reads
READS=`awk -F, '$1 == "illumina-18_S15" {print $2}' ./results/removal_summary.csv`
if [[ "$READS" != "2" ]]; then 
    echo "Incorrect output: Number of Human Reads"
    echo "  Expected: 2, Got: $READS"
    exit 1
fi

# 2. Total Reads Kept
READS=`awk -F, '$1 == "illumina-4_S12" {print $4}' ./results/removal_summary.csv`
if [[ "$READS" != "100" ]]; then 
    echo "Incorrect output: Number of Paired Reads Kept"
    echo "  Expected: 100, Got: $READS"
    exit 1
fi

# 3. Downsample Maximum
READS=`awk -F, '$1 == "illumina-4_S12" {print $6}' ./results/removal_summary.csv`
if [[ "$READS" != "100" ]]; then 
    echo "Incorrect output: Downsampling Maximum"
    echo "  Expected: 100, Got: $READS"
    exit 1
fi

# 4. Downsample Seed
SEED=`awk -F, '$1 == "illumina-4_S12" {print $7}' ./results/removal_summary.csv`
if [[ "$SEED" != "42" ]]; then 
    echo "Incorrect output: Downsampling Seed"
    echo "  Expected: 42, Got: $SEED"
    exit 1
fi

# 5. Region count file exists
if [[ ! -f "./results/downsample_stats/illumina-4_S12_region_counts.csv" ]]; then
    echo "illumina-4_S12_region_counts.csv expected output file not found"
    exit 1
fi

# Reset and Track
mv .nextflow.log artifacts/illumina_downsampled.nextflow.log
rm -rf results work/ .nextflow*

echo "Done"
