#!/usr/bin/env python
# encoding: utf-8

"""
screen_vcf_files.py

Created by Brant Faircloth on 16 December 2013.
Copyright 2013 Brant C. Faircloth. All rights reserved.

"""

import os
import vcf
import argparse

import pdb


def get_args():
    parser = argparse.ArgumentParser(description='Parse fastq files and drop reads containing Ns.')
    parser.add_argument(
            '--input',
            required=True,
            help="The vcf file to process"
        )
    parser.add_argument(
            '--output',
            required=True,
            help="The vcf file to write"
        )
    parser.add_argument(
            '--complete',
            type=float,
            default=100.00,
            help="The percentage of missing data to allow."
        )
    return parser.parse_args()

def main():
    args = get_args()
    count = 0
    with open(args.input, 'r') as infile:
        with open(args.output, 'w') as outfile:
            vcf_reader = vcf.Reader(infile)
            vcf_writer = vcf.Writer(outfile, vcf_reader)
            for record in vcf_reader:
                complete = []
                # make sure that we get only passing (just in case)
                if record.FILTER == []:
                    for sample in record.samples:
                        # check for min completeness
                        if sample.data.GT == None:
                            complete.append(0)
                        else:
                            complete.append(1)
                    if (100 * (sum(complete) / float(len(complete))) >= args.complete):
                        vcf_writer.write_record(record)
                print count
                count += 1

if __name__ == '__main__':
    main()
