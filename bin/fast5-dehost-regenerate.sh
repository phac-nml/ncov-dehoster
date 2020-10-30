#!/bin/bash

# Set up to re-generate dehosted fast5 files after the re-basecalling step
# Done in this bash script as NF does not like the awk line

# Barcode folder of dehosted fastq files
dehosted_fastq_barcode=$1

# Barcode folder name (ex. barcode21 or unclassified)
barcodeName=$2

# Fast5_dehosted combined directory from nanostripper with all of the original dehosted fast5 files
fast5_dehosted_only=$3

# Get the fastq headers
awk -F ' ' '(NR%4==1) {print $1}' ${dehosted_fastq_barcode}/* | sed 's/@//g'  > ${barcodeName}.txt

# Make fast5 dehosted directory
mkdir -p fast5_pass/$barcodeName

# Run ONT fast5 subset to make the new fast5 files
fast5_subset  -t 16 -i $fast5_dehosted_only -f $barcodeName'_' -s fast5_pass/$barcodeName -l ${barcodeName}.txt --recursive
