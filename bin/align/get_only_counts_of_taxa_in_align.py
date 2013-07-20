#!/usr/bin/env python
# encoding: utf-8

"""
get_counts_of_taxa_in_align.py

Created by Brant Faircloth on 16 August 2011.
Copyright 2011 Brant C. Faircloth. All rights reserved.

The program iterates through a folder of nexus files and returns the count
of alignments having more than "--percent" of "--taxa" taxa.
"""


import os
import glob
import math
import shutil
import argparse
from Bio import AlignIO
from phyluce.helpers import is_dir

#import pdb


def get_args():
    parser = argparse.ArgumentParser(description='Match UCE probes to assembled contigs and store the data')
    parser.add_argument('nexus', help='The directory containing the nexus files', type=is_dir)
    parser.add_argument('taxa', type=int, help='The number of taxa expected')
    parser.add_argument('output', type=is_dir, help='The output dir in which to store copies of the alignments')
    parser.add_argument('--percent', type=float, default=0.5, help='The percent of taxa to require')
    return parser.parse_args()


def get_files(input_dir):
    return glob.glob(os.path.join(os.path.expanduser(input_dir), '*.nex*'))


def main():
    args = get_args()
    # iterate through all the files to determine the longest alignment
    files = get_files(args.nexus)
    min_count = int(math.floor(args.percent * args.taxa))
    counts = 0
    for f in files:
        aln = AlignIO.read(f, 'nexus')
        if len(aln) < min_count:
            pass
        else:
            counts += 1
            shutil.copyfile(f, os.path.join(args.output, os.path.basename(f)))
    print "Copied {0} alignments containing ≥ {1} proportion of taxa (n = {2})".format(
        counts,
        args.percent,
        min_count
    )

if __name__ == '__main__':
    main()
