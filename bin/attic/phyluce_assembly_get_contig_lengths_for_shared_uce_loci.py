#!/usr/bin/env python
# encoding: utf-8

"""
get_contig_lengths_for_uce_loci.py

Created by Brant Faircloth on 02 July 2011.
Copyright 2011 Brant C. Faircloth. All rights reserved.
"""

import os
import sys
import argparse
from collections import defaultdict
from seqtools.sequence import fasta

import pdb

def get_args():
    parser = argparse.ArgumentParser(description='Get the average contig length, by species, within a combined fasta')
    parser.add_argument('fasta', help='The fasta input')
    parser.add_argument('--output', help = 'The output file',
        type=argparse.FileType('w'), default=sys.stdout)
    return parser.parse_args()

def main():
    args = get_args()
    records = fasta.FastaReader(args.fasta)
    lengths = defaultdict(list)
    for sequence in records:
        # BEWARE:  this may cause name clash, which will error out
        org = sequence.identifier.split(' ')[0].split('_')[-2]
        lengths[org].append(len(sequence))
    for org,l in lengths.iteritems():
        #pdb.set_trace()
        args.output.write("{0}\t{1}\n".format(org, float(sum(l))/len(l)))

if __name__ == '__main__':
    main()
