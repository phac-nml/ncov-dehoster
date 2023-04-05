#!/bin/bash
set -eo pipefail

# Create a Cache Dir
mkdir -p conda_cache_dir

# --------------------------------------------------------------------- #
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
    echo "Incorrect output: Number of Poor Reads"
    echo "  Expected: 2, Got: $READS"
    exit 1
fi

# 2. Total Reads Kept
READS=`awk -F, '$1 == "nanopore-2" {print $4}' ./results/test-1-minimap2-flat/removal_summary.csv`
if [[ "$READS" != "2013" ]]; then 
    echo "Incorrect output: Number of Reads Kept"
    echo "  Expected: 2013, Got: $READS"
    exit 1
fi

# 3. Check Formatting
FILE_COUNT=`ls ./results/test-1-minimap2-flat/run/fastq_pass/ | grep fastq -c`
if [[ "$FILE_COUNT" != "2" ]]; then 
    echo "Incorrect output: Number of fastq files does not appear to be correct"
    echo "  Expected: 2, Got: $FILE_COUNT"
    exit 1
fi

# Reset and Track
mv .nextflow.log artifacts/minimap2_flat.nextflow.log
rm -rf results work/ .nextflow*

# --------------------------------------------------------------------- #
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
    echo "Incorrect output: Number of Poor Reads"
    echo "  Expected: 2, Got: $READS"
    exit 1
fi

# 2. Total Reads Kept
READS=`awk -F, '$1 == "nanopore-2" {print $4}' ./results/test-2-minimap2-partial-human/removal_summary.csv`
if [[ "$READS" != "2013" ]]; then 
    echo "Incorrect output: Number of Reads Kept"
    echo "  Expected: 2013, Got: $READS"
    exit 1
fi

# 3. Check Formatting
DIR_COUNT=`ls -1d ./results/test-2-minimap2-partial-human/run/fastq_pass/*/ | wc -l`
if [[ "$DIR_COUNT" != "2" ]]; then 
    echo "Incorrect output: Number of fastq directories does not appear to be correct"
    echo "  Expected: 2, Got: $DIR_COUNT"
    exit 1
fi

# Reset and Track
mv .nextflow.log artifacts/minimap2_expanded.nextflow.log
rm -rf results work/ .nextflow*

# --------------------------------------------------------------------- #
### Run flat and test for downsampling ###
nextflow run ./main.nf \
    -profile conda,test \
    --cache ./conda_cache_dir \
    --nanopore \
    --minimap2 \
    --downsample \
    --downsample_count 100 \
    --downsample_seed 42 \
    --fastq_directory $PWD/.github/data/nanopore/ \
    --run_name 'test-3-minimap2-partial-human-flat-down100' \
    --flat \
    --human_ref $PWD/.github/data/partial_hg38_ref.fa

### Check Outputs ###
# 1. Num Poor Reads
READS=`awk -F, '$1 == "nanopore-1" {print $3}' ./results/test-3-minimap2-partial-human-flat-down100/removal_summary.csv`
if [[ "$READS" != "2" ]]; then 
    echo "Incorrect output: Number of Poor Reads"
    echo "  Expected: 2, Got: $READS"
    exit 1
fi

# 2. Total Reads Kept
READS=`awk -F, '$1 == "nanopore-2" {print $4}' ./results/test-3-minimap2-partial-human-flat-down100/removal_summary.csv`
if [[ "$READS" != "100" ]]; then 
    echo "Incorrect output: Number of Reads Kept"
    echo "  Expected: 100, Got: $READS"
    exit 1
fi

# 3. Check Formatting
FILE_COUNT=`ls ./results/test-3-minimap2-partial-human-flat-down100/run/fastq_pass/ | grep fastq -c`
if [[ "$FILE_COUNT" != "2" ]]; then 
    echo "Incorrect output: Number of fastq files does not appear to be correct"
    echo "  Expected: 2, Got: $FILE_COUNT"
    exit 1
fi

# 4. Check Downsampling Max
READS=`awk -F, '$1 == "nanopore-1" {print $6}' ./results/test-3-minimap2-partial-human-flat-down100/removal_summary.csv`
if [[ "$READS" != "100" ]]; then 
    echo "Incorrect output: Downsampling Maximum"
    echo "  Expected: 100, Got: $READS"
    exit 1
fi

# 5. Check Downsampling Seed
SEED=`awk -F, '$1 == "nanopore-1" {print $7}' ./results/test-3-minimap2-partial-human-flat-down100/removal_summary.csv`
if [[ "$SEED" != "42" ]]; then 
    echo "Incorrect output: Downsampling Seed"
    echo "  Expected: 42, Got: $SEED"
    exit 1
fi

# Reset and Track
mv .nextflow.log artifacts/minimap2_expanded_flat_downsample.nextflow.log
rm -rf results work/ .nextflow*

# --------------------------------------------------------------------- #
### Run non-flat and test for downsampling ###
nextflow run ./main.nf \
    -profile conda,test \
    --cache ./conda_cache_dir \
    --nanopore \
    --minimap2 \
    --downsample \
    --downsample_count 100 \
    --downsample_seed 42 \
    --fastq_directory $PWD/.github/data/nanopore/ \
    --run_name 'test-4-minimap2-partial-human-down100' \
    --human_ref $PWD/.github/data/partial_hg38_ref.fa

### Check Outputs ###
# 1. Num Poor Reads
READS=`awk -F, '$1 == "nanopore-1" {print $3}' ./results/test-4-minimap2-partial-human-down100/removal_summary.csv`
if [[ "$READS" != "2" ]]; then 
    echo "Incorrect output: Number of Poor Reads"
    echo "  Expected: 2, Got: $READS"
    exit 1
fi

# 2. Total Reads Kept
READS=`awk -F, '$1 == "nanopore-2" {print $4}' ./results/test-4-minimap2-partial-human-down100/removal_summary.csv`
if [[ "$READS" != "100" ]]; then 
    echo "Incorrect output: Number of Reads Kept"
    echo "  Expected: 100, Got: $READS"
    exit 1
fi

# 3. Check Formatting
DIR_COUNT=`ls -1d ./results/test-4-minimap2-partial-human-down100/run/fastq_pass/*/ | wc -l`
if [[ "$DIR_COUNT" != "2" ]]; then 
    echo "Incorrect output: Number of fastq directories does not appear to be correct"
    echo "  Expected: 2, Got: $DIR_COUNT"
    exit 1
fi

# 4. Check Downsampling Max
READS=`awk -F, '$1 == "nanopore-1" {print $6}' ./results/test-4-minimap2-partial-human-down100/removal_summary.csv`
if [[ "$READS" != "100" ]]; then 
    echo "Incorrect output: Downsampling Maximum"
    echo "  Expected: 100, Got: $READS"
    exit 1
fi

# 5. Check Downsampling Seed
SEED=`awk -F, '$1 == "nanopore-1" {print $7}' ./results/test-4-minimap2-partial-human-down100/removal_summary.csv`
if [[ "$SEED" != "42" ]]; then 
    echo "Incorrect output: Downsampling Seed"
    echo "  Expected: 42, Got: $SEED"
    exit 1
fi

# Reset and Track
mv .nextflow.log artifacts/minimap2_expanded_downsample.nextflow.log
rm -rf results work/ .nextflow*

echo "Done"
