#!/usr/bin/env python
# encoding: utf-8
"""
File: get_fastq_lengths.py
Author: Brant Faircloth

Created by Brant Faircloth on 21 July 2012 12:07 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: 

"""

import gzip
import argparse


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Get summary (length) data from fastq""")
    parser.add_argument(
            "fastq",
            help="""The fastq file to summarize"""
        )
    parser.add_argument(
            "--csv",
            action="store_true",
            default=False,
            help="""Give output in CSV"""
        )
    return parser.parse_args()


def main():
    args = get_args()
    if args.fastq.endswith('.gz'):
        f = gzip.open(args.fastq)
    else:
        f = open(args.fastq, 'rU')
    count = []
    while True:
        l1, l2, l3, l4 = [f.readline() for line in range(4)]
        if not l2:
            break
        count.append(len(l2))
    if not args.csv:
        print "Reads:\t\t{:,}".format(len(count))
        print "Bp:\t\t{:,}".format(sum(count))
        print "Avg. len:\t{:,}".format(sum(count) / len(count))
        print "Min. len:\t{:,}".format(min(count))
        print "Max. len:\t{:,}".format(max(count))
    else:
        print "{},{},{},{}".format(args.fastq, len(count), sum(count), sum(count) / len(count))



if __name__ == '__main__':
    main()