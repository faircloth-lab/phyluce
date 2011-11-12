#!/usr/bin/env python
# encoding: utf-8

"""
find_snps_in_bed_interval.py

Created by Brant Faircloth on 11 June 2011.
Copyright 2011 Brant C. Faircloth. All rights reserved.
"""

import os
import sys
import argparse
from operator import itemgetter

import pdb

def get_args():
    parser = argparse.ArgumentParser(description='Match UCE probes to assembled contigs and store the data')
    parser.add_argument('uce', help='The UCE BED file', type=argparse.FileType('rU'))
    parser.add_argument('snp', help='The SNP intersection BED file', type=argparse.FileType('rU'))
    parser.add_argument('--output', help = 'The output file',
        type=argparse.FileType('w'), default=sys.stdout)
    return parser.parse_args()
    
def get_intersect(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))
    
def fix_ints(line, splitchar = '\t'):
    lsplit = line.strip('\n').split(splitchar)
    return [lsplit[0], int(lsplit[1]), int(lsplit[2]), lsplit[3]]

def main():
    args = get_args()
    snps = [fix_ints(line) for line in args.snp if not line.startswith('track')]
    mapped = []
    for line in args.uce:
        if not line.startswith('track'):
            chromo, start, end, uce = fix_ints(line, ' ')
            for pos, snp in enumerate(snps):
                if snp[0] == chromo:
                    if get_intersect([start, end], [snp[1], snp[2]]):
                        mapped.append([uce, chromo, start, end, snp[3], snp[1], snp[2]])
    args.output.write('UCE,chromo,uce-start,uce-end,snp-name,snp-start,snp-end\n')
    for i in sorted(mapped, key = itemgetter(0)):
        args.output.write("{},{},{},{},{},{},{}\n".format(i[0], i[1], i[2], i[3], i[4], i[5], i[6]))

if __name__ == '__main__':
    main()