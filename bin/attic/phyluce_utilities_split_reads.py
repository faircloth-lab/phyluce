#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2015 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 01 March 2012 22:03 PST (-0800)
"""

import sys
import argparse
from seqtools.sequence import fastq

def get_args():
    parser = argparse.ArgumentParser(description='Split an interleaved, paired-end fastq file into two files')
    parser.add_argument(
        '--input',
        required=True,
        help='The input fastq file'
    )
    parser.add_argument(
        '--read1',
        required=True,
        help='The output read1 fastq file name'
    )
    parser.add_argument(
        '--read2',
        required=True,
        help='The output read2 fastq file name'
    )
    args = parser.parse_args()
    return args

def main():
    args = get_args()
    r1 = fastq.FasterFastqWriter(args.read1)
    r2 = fastq.FasterFastqWriter(args.read2)
    # read all of our files into fastq iterators
    rc = 0
    reads = fastq.FasterFastqReader(args.input)
    sys.stdout.write("Splitting reads (1 dot = 10,000 pairs): ")
    for read in reads:
        if read[0].split(' ')[1].split(':')[0] == '1':
            r1.write(read)
            first = read[0].split(' ')[0]
        else:
            assert first == read[0].split(' ')[0], "File does not appear interleaved."
            r2.write(read)
            if rc != 0 and rc % 10000 == 0:
                sys.stdout.write('.')
                sys.stdout.flush()
            rc += 1
    print ""
    reads.close()
    r1.close()
    r2.close()


if __name__ == '__main__':
    main()
