#!/usr/bin/env python3

import argparse
import pysam

def init_parser():
    '''
    Parser Arguments to pass to script from CL
    '''

    parser = argparse.ArgumentParser()
    parser.add_argument('-k', '--keep_sam', required=True, help='SAM file of mapped reads to keep')
    parser.add_argument('-r', '--remove_sam', required=True, help='SAM file of mapped reads to remove')
    parser.add_argument('-q', '--keep_minimum_quality', required=False, type=int, default=60, help='Minimum quality of the reads to keep')
    parser.add_argument('-Q', '--remove_minimum_quality', required=False, type=int, default=0, help='Minimum quality of the reads to be included in removal')
    parser.add_argument('-m', '--minimum_read_length', required=False, type=int, default=251, help='Minimum length of reads to keep')
    parser.add_argument('-M', '--maximum_read_length', required=False, type=int, default=0, help='Maximum length of reads to keep')
    parser.add_argument('-o', '--output', required=False, default='out.sam', help='Output sam name')

    return parser


def get_reads_to_remove(samfile_path, input_mapping_quality, reads_to_remove=set(), count=0):

    samfile = pysam.AlignmentFile(samfile_path, "r")

    for read in samfile.fetch():
        if read.mapping_quality >= input_mapping_quality:
            reads_to_remove.add(read.query_name)
            count += 1
    samfile.close()

    return reads_to_remove, count


def generate_dehosted_sam(samfile_in, samfile_out, input_mapping_quality, removal_read_ids, reads_found=[]):

    samfile = pysam.AlignmentFile(samfile_in, "r")
    header = samfile.header

    for read in samfile.fetch():
        if read.query_name in removal_read_ids:
            pass

        elif read.mapping_quality >= input_mapping_quality:
            reads_found.append(read)

    samfile.close()

    with pysam.AlignmentFile(samfile_out, "w", header=header) as outfile:
        for read in reads_found:
            outfile.write(read)


def main():
    parser = init_parser()
    args = parser.parse_args()

    remove_reads_set, num_removed = get_reads_to_remove(args.remove_sam, args.remove_minimum_quality)

    generate_dehosted_sam(args.keep_sam, args.output, args.keep_minimum_quality, remove_reads_set)


if __name__ == "__main__":
    main()
