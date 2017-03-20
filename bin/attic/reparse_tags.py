#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2017 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 17 March 2017 15:15 CDT (-0500)
"""

import pdb
import sys
import gzip
import argparse
from Bio import SeqIO

from collections import Counter


def get_args():
    parser = argparse.ArgumentParser(description='Parse fastqs from input based on tag sequence')
    parser.add_argument('--input', nargs='?', default=sys.stdin)
    parser.add_argument('--tag', required = True, dest = 'tag')
    parser.add_argument('--output', nargs='?', default=sys.stdout)
    return parser.parse_args()

def main():
    args = get_args()
    count = 0
    with gzip.open(args.output, 'wb') as outfile:
        with gzip.open(args.input, 'rb') as infile:
            fastqs = SeqIO.parse(infile, 'fastq')
            for read in fastqs:
                read_tag = read.description.split(' ')[1].split(":")[-1]
                if read_tag == args.tag:
                    outfile.write(read.format('fastq'))
                count += 1
                if count%1000000 == 0:
                    print count
                #pdb.set_trace()

if __name__ == '__main__':
    main()
