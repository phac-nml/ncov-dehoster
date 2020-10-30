#!/bin/bash

# Set the numbr of parallel runs to do for basecalling
parallel_basecalling=$1

# Sets the number of gpus on the node
# Currently not working, have to find a way to put into cuda part
gpus=$2

for folder in `ls -1 -d */`; do echo $folder; done | parallel -j $parallel_basecalling   guppy_basecaller -c dna_r9.4.1_450bps_hac.cfg -r -i {} -s fastq_pass_dehosted_only/{} -x \"cuda:'{= $_=$job->slot()%2=}'\"
