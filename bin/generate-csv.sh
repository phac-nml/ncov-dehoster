#!/bin/bash

# Setting names nand paths to generate the csv file
# This is a bash script so we can use bash variables easier

# Need barcode name to keep track of the files
barcodeName=$1

# Path to the nanostripper summary
nanostripper_summary_path=$2

# Path the the newly made filename_mapping file
filename_mapping_path=$3

# Pipeline revision has to keep track of what the dehosting script looked like in case of issue
rev=$4

# Set up csv file with headers
echo "barcode,num_stripped_reads,num_kept_reads,pipeline_revision" > ${barcodeName}.csv

# Get values
stripped=`wc -l < ${nanostripper_summary_path}`
kept=`wc -l < ${filename_mapping_path}`

# Generate file
echo "${barcodeName},${stripped},${kept},${rev}" >> ${barcodeName}.csv
