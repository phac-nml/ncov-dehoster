#!/bin/bash

#
# Downsample amplicon sequencing BAM file with samtools #
#

### Functions ###
# Check if file is valid
valid_file() {
    if [ ! -f $1 ]; then
        echo "ERROR: Given arg '$2 $1' is not a valid file"
        exit 1
    fi
}
# Check if int is valid 
valid_int() {
    if [[ $1 != +([0-9]) ]]; then
        echo "ERROR: Given arg '$2 $1' is not an integer!"
        exit 1
    fi
}
# Check if platform is available
valid_platform() {
    if [ "$1" == "illumina" -o "$1" == "nanopore" ]; then
        :
    else 
        echo "ERROR: Given '$2 $1' is not a valid platform. Please pass either 'illumina' or 'nanopore'"
        exit 1
    fi
}

# Help
HELP="
USAGE: ./downsample_bam.sh -b BAM -a BED -p [illumina nanopore]

REQUIRED:
    -b  --bam        : Path - input BAM file to downsample
    -a  --amplicons  : Path - amplicon bed file to
    -p  --platform   : String - Platform used. Options are 'illumina' and 'nanopore'
    
OPTIONAL:
    -c  --read-count : Int - Target output read count per amplicon. Default: 1200
    -s  --seed       : Int - Samtools seed to downsample with. Default: 101
    --sample-name    : String - Sample name for outputs files. Default: Filename cut at first '.'
"

### Inputs ###
##############
# Required files and default variables
BAM_IN=""
AMPLICON_IN=""
PLATFORM=""
DOWNSAMPLE_SEED=101
TARGET_AMPLICON_COUNT=1200
SAMPLE_NAME=""
PAIRED_FLAG=""

# Check for Args #
if [ $# -eq 0 ]; then
    echo "$HELP"
    exit 0
fi

# Set and Check Arguments #
while [ $# -gt 0 ];
do
    if [ $1 = "--bam" -o $1 = "-b" ]; then
        shift
        BAM_IN=$1
        valid_file $BAM_IN '--bam'
        shift
    elif [ $1 = "--amplicons" -o $1 = "-a" ]; then
        shift
        AMPLICON_IN=$1
        valid_file $AMPLICON_IN '--amplicons'
        shift
    elif [ $1 = "--platform" -o $1 = "-p" ]; then
        shift
        PLATFORM=$1
        valid_platform $PLATFORM '--platform'
        shift
    elif [ $1 = "--read-count" -o $1 = "-c" ]; then
        shift
        TARGET_AMPLICON_COUNT=$1
        valid_int $TARGET_AMPLICON_COUNT '--read-count'
        shift
    elif [ $1 = "--seed" -o $1 = "-s" ]; then
        shift
        DOWNSAMPLE_SEED=$1
        valid_int $DOWNSAMPLE_SEED '--seed'
        shift
    elif [ $1 = "--sample-name" ]; then
        shift
        SAMPLE_NAME=$1
        shift
    else
        echo "ERROR: $1 is not a known argument"
        exit 1
    fi
done

# Check for required inputs
if [[ "$BAM_IN" == "" || "$AMPLICON_IN" == "" ]]; then
    echo "Missing input file"
    exit 1
fi

# Set default values
if [[ "$SAMPLE_NAME" == "" ]]; then
    SAMPLE_NAME=`echo $BAM_IN | cut -f1 -d'.'`
fi
if [[ "$PLATFORM" == "illumina" ]]; then
    PAIRED_FLAG="-f3"
fi

# Change target amplicon count to an even number
if (( TARGET_AMPLICON_COUNT % 2 == 1 )); then
    ((TARGET_AMPLICON_COUNT++))
fi

### Script ###
##############
# Create regions file to use to downsample.
#   Input bed file first 3 columns `contigID	ampliconStart	ampliconEnd` and the rest are whatever
#   Format is a samtools region string followed by the start and end of the amplicon
awk -F"\t" 'NR >0 {print $1":" $2"-"$3","$2","$3}' $AMPLICON_IN > regions.txt

# Setup output tracking
echo samtools_region,initial_read_count,expected_final_read_count,fraction_kept,seed > ${SAMPLE_NAME}_region_counts.csv

# Header for samtools later steps
samtools view -H ${BAM_IN} > header_line.sam

# For each region, we want to get all reads that are at least 60% in the amplicon
#   This is to not count any reads twice with the tiling
while read region_info; do
    samtools_region=`echo ${region_info} | cut -f1 -d','`
    region_min=`echo ${region_info} | cut -f2 -d','`
    region_max=`echo ${region_info} | cut -f3 -d','`
    echo "$samtools_region"

    # Filter for reads in the amplicon region with awk and based on input platform. Threshold is >= 60% in amplicon region
    #   Illumina: Check if read is forward or reverse based on the total region lengths sign and then calculate the boundries based on that
    #   Nanopore: Calculate the boundries for the read based simply on the start position and read length
    samtools view $PAIRED_FLAG $BAM_IN $samtools_region \
        | awk -v region_min="$region_min" -v region_max="$region_max" -v platform="$PLATFORM" '{
            threshold=0.60 ;
            if(platform=="illumina"){
                total_len=$9 ;
                if(total_len > 0){
                    low_bound=$4 ;
                    up_bound=$4+$9 ;
                }
                else {
                    low_bound=$8 ;
                    up_bound=$8+(-1*$9) ;
                    total_len=(-1*$9) ;
                }
            }
            else {
                total_len=length($10)
                low_bound=$4 ;
                up_bound=$4+total_len ;
            }
            (up_bound < region_max) ? min_up=up_bound : min_up=region_max ;
            (low_bound > region_min) ? max_low=low_bound : max_low=region_min ;
            overlap=min_up-max_low ;
            frac=overlap/total_len ;
            if( frac >= threshold ) print $0
        }' > ${region_min}-${region_max}.tmp.sam

    # Check if number of reads is less than amplicon target
    count=`wc -l <${region_min}-${region_max}.tmp.sam`
    if [ "$count" -lt $TARGET_AMPLICON_COUNT ]; then
        cat header_line.sam $region_min-$region_max.tmp.sam \
            | samtools view $PAIRED_FLAG -bh \
            > ${region_min}-${region_max}_downsampled.tmp.bam
        echo $samtools_region,$count,$count,1.00,$DOWNSAMPLE_SEED >> ${SAMPLE_NAME}_region_counts.csv
        continue
    fi

    # Samtools utilizes `-s seed.fraction` for its subsampling so want a fair number of digits to get close to input wanted read count
    frac=`awk -v seed=$DOWNSAMPLE_SEED -v want=$TARGET_AMPLICON_COUNT -v total=$count 'BEGIN {printf("%.0f",(want/total)*10000)}'`
    cat header_line.sam $region_min-$region_max.tmp.sam \
        | samtools view $PAIRED_FLAG -bhs ${seed}.${frac} \
        > ${region_min}-${region_max}_downsampled.tmp.bam

    # Track - Note that the final read count will not be found at each position in the amplicon as some reads might only contain a part of it
    frac="0.$frac"
    new_count=`awk -v frac=$frac -v total=$count 'BEGIN {printf("%.0f",frac*total)}'`
    echo $samtools_region,$count,$new_count,$frac,$DOWNSAMPLE_SEED >> ${SAMPLE_NAME}_region_counts.csv

done < regions.txt

# Merge
samtools merge ${SAMPLE_NAME}_downsampled.bam *_downsampled.tmp.bam

# Remove unneeded files
rm *.tmp.sam
rm header_line.sam
rm *_downsampled.tmp.bam
echo "done"
