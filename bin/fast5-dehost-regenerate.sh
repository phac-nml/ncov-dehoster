#!/bin/bash

# Set up to re-generate dehosted fast5 files after the re-basecalling step
# Done in this bash script as NF does not like the awk line

### Inputs ###
# Barcode folder of dehosted fastq files
DEHOSTED_BARCODE_FASTQ_FOLDER=$1
# Barcode folder name (ex. barcode21 or unclassified)
BARCODE_NAME=$2
# Fast5_dehosted combined directory from nanostripper with all of the original dehosted fast5 files
FAST5_IN=$3
# Threads to use for fast5_subset
THREADS=$4

### Run ###
# Get the fastq headers
awk -F ' ' '(NR%4==1) {print $1}' ${DEHOSTED_BARCODE_FASTQ_FOLDER}/* | sed 's/@//g'  > ${BARCODE_NAME}.txt

# Make fast5 dehosted directory for the barcode
mkdir -p fast5_pass/$BARCODE_NAME

# Run ONT fast5 subset to make the new fast5 files
fast5_subset  -t $THREADS -i $FAST5_IN -f $BARCODE_NAME'_' -s fast5_pass/$BARCODE_NAME -l ${BARCODE_NAME}.txt --recursive
