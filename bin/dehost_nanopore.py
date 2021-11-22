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
    parser.add_argument('-q', '--keep_minimum_quality', required=False, type=int, default=60, help='Minimum quality of the reads to keep. Default: 60')
    parser.add_argument('-Q', '--remove_minimum_quality', required=False, type=int, default=10, help='Minimum quality of the reads to be included in removal. Default: 10')
    parser.add_argument('-o', '--output', required=False, default='out.bam', help='Output BAM name')

    return parser

def get_reads_to_remove(bamfile_path, contig_ID, remove_minimum_quality, keep_minimum_quality, h_count=0, p_count=0, reads_to_remove_set=set(), read_list=[]):
    # File input and run through
    bamfile = pysam.AlignmentFile(bamfile_path, "rb")
    header = bamfile.header
    # Removal of host reads if they are above threshold
    for read in bamfile.fetch():
        if read.reference_name != contig_ID and read.mapping_quality >= remove_minimum_quality:
            reads_to_remove_set.add(read.query_name)
            h_count += 1

    # Have to go through again as we need to remove reads that have mapped twice, can probably do it in one pass but blanking on a good solution
    for read in bamfile.fetch():
        # Already found human reads, don't add them
        if read.query_name in reads_to_remove_set:
            pass
        # Keep reads that are above the input minimum mapping quality and not human
        elif read.mapping_quality >= keep_minimum_quality:
            read_list.append(read)
        # Poor quality
        else:
            p_count += 1
    bamfile.close()

    return read_list, h_count, p_count, header


def generate_dehosted_output(read_list, keep_header, bamfile_out):
    with pysam.AlignmentFile(bamfile_out, "w", header=keep_header) as outfile:
        for read in read_list:
            outfile.write(read)

def main():
    parser = init_parser()
    args = parser.parse_args()

    sample_name = os.path.splitext(Path(args.file).stem)[0]

    keep_read_list, h_count, p_count, header = get_reads_to_remove(args.file, args.keep_id, args.remove_minimum_quality, args.keep_minimum_quality)

    generate_dehosted_output(keep_read_list, header, args.output)

    if len(keep_read_list) == 0:
        percentage_kept = 0
    else:
        percentage_kept = len(keep_read_list)/(h_count + p_count + len(keep_read_list)) * 100

    line = {    'sample' : sample_name,
                'human_reads_filtered' : h_count, 
                'poor_quality_reads_filtered' : p_count,
                'paired_reads_kept' : len(keep_read_list),
                'percentage_kept' : "{:.2f}".format(percentage_kept)
            }

    with open('{}_stats.csv'.format(sample_name), 'w') as csvfile:
        header = line.keys()
        writer = csv.DictWriter(csvfile, fieldnames=header)
        writer.writeheader()
        writer.writerow(line)

if __name__ == "__main__":
    main()
