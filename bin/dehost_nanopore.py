#!/usr/bin/env python3

import argparse
import csv
import pysam
import os

from pathlib import Path

def init_parser():
    '''
    Parser Arguments to pass to script from CL
    '''

    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--file', required=True, help='Composite Reference BAM file of mapped reads')
    parser.add_argument('-k', '--keep_id', required=False, default='MN908947.3', type=str, help='Reference ID of genome to keep. Default: MN908947.3')
    parser.add_argument('-m', '--min_reads', required=False, default=1, type=int, help='Minimum number of reads required to generate dehosted BAM file')
    parser.add_argument('-q', '--keep_minimum_quality', required=False, type=int, default=60, help='Minimum quality of the reads to keep. Default: 60')
    parser.add_argument('-Q', '--remove_minimum_quality', required=False, type=int, default=10, help='Minimum quality of the reads to be included in removal. Default: 10')
    parser.add_argument('-o', '--output', required=False, default='out.bam', help='Output BAM name')
    parser.add_argument('-R', '--revision', required=False, default='NA', help='Pass a pipeline commit hash to keep track of what version was ran')
    parser.add_argument('--downsampled', required=False, action='store_true', help='Pass if data is downsampled')
    parser.add_argument('--downsampled_count', required=False, type=int, default=100000, help='Maximum count of reads to keep when downsampling')
    parser.add_argument('--downsampled_seed', required=False, type=int, default=100, help='Downsample seed used')

    return parser

def keep_reads_by_contig_id(bamfile_path, contig_ID, remove_minimum_quality, keep_minimum_quality, h_count=0, p_count=0, reads_to_remove_set=set(), read_list=[]):
    '''
    PURPOSE:
        Keep reads that are above the minimum quality threshold and match the wanted contig_ID.
            Tracks the count of human and poor quality reads
            Saves the BAM header for output later
    INPUTS:
        - bamfile_path: Path to the input BAM file to dehost
        - contig_ID: String of the contig ID value to keep
        - remove_minimum_quality: Int value of the minimum read value to remove
        - keep_minimum_quality: Int value of the minimum read value to retain
    RETURNS:
        - read_list: List of pysam AlignedSegments that match wanted qualities
        - h_count: Int count of human reads removed
        - p_count: Int count of poor quality reads removed
        - header: Pysam AlignmentHeader
    '''
    # File input and save the BAM header
    bamfile = pysam.AlignmentFile(bamfile_path, "rb")
    header = bamfile.header

    # Selection of reads based on mapping quality and what contig they map to
    for read in bamfile.fetch():
        # Don't keep secondary reads, no need to capture them either
        if read.is_secondary:
            pass
        # Remove reads that don't map to our expected reference ID
        elif read.reference_name != contig_ID and read.mapping_quality >= remove_minimum_quality:
            reads_to_remove_set.add(read.query_name)
            h_count += 1
        # Keep reads that are above the input minimum mapping quality and not already a human read (sort has human mapping reads first)
        elif read.mapping_quality >= keep_minimum_quality and read.query_name not in reads_to_remove_set:
            read_list.append(read)
        # Everything else is poor quality and removed
        else:
            p_count += 1
    bamfile.close()

    return read_list, h_count, p_count, header


def generate_dehosted_output(read_list, keep_header, bamfile_out):
    '''
    PURPOSE:
        Generate output dehosted BAM file to output path
    INPUTS:
        - read_list: List of pysam AlignedSegments that match wanted qualities
        - keep_header: Pysam AlignmentHeader from input BAM file
        - bamfile_out: String path to output file
    '''
    with pysam.AlignmentFile(bamfile_out, "w", header=keep_header) as outfile:
        for read in read_list:
            outfile.write(read)

def main():
    # Get args set to go
    parser = init_parser()
    args = parser.parse_args()

    # Capture the sample name from the file
    sample_name = os.path.splitext(Path(args.file).stem)[0]

    # Keep reads based on input contig ID and mapping qualities given and then generate the output bam file using input files header with pysam
    keep_read_list, h_count, p_count, header = keep_reads_by_contig_id(args.file, args.keep_id, args.remove_minimum_quality, args.keep_minimum_quality)
    kept_count = len(keep_read_list)
    if kept_count >= args.min_reads:
        generate_dehosted_output(keep_read_list, header, args.output)
        output_generated = True
    else:
        output_generated = False

    # Set up the output CSV file for tracking
    if kept_count == 0:
        percentage_kept = 0
    else:
        percentage_kept = kept_count/(h_count + p_count + kept_count) * 100

    # Output based on if downsampled or not
    if args.downsampled:
        if kept_count > args.downsampled_count:
            kept_count = args.downsampled_count
            percentage_kept = kept_count/(h_count + p_count + kept_count) * 100

        line = {    'sample' : sample_name,
                    'human_reads_filtered' : h_count, 
                    'poor_quality_reads_filtered' : p_count,
                    'reads_kept' : kept_count,
                    'percentage_kept' : "{:.2f}".format(percentage_kept),
                    'downsample_maximum_reads' : args.downsampled_count,
                    'downsample_seed' : args.downsampled_seed,
                    'github_commit' : args.revision}
    else:
        line = {    'sample' : sample_name,
                    'human_reads_filtered' : h_count, 
                    'poor_quality_reads_filtered' : p_count,
                    'reads_kept' : kept_count,
                    'percentage_kept' : "{:.2f}".format(percentage_kept),
                    'meets_count_filter' : output_generated,
                    'github_commit' : args.revision
                }

    # Write out the output CSV file
    with open('{}_stats.csv'.format(sample_name), 'w') as csvfile:
        header = line.keys()
        writer = csv.DictWriter(csvfile, fieldnames=header)
        writer.writeheader()
        writer.writerow(line)


if __name__ == "__main__":
    main()
