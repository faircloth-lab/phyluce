#!/usr/bin/env python
# encoding: utf-8

"""
get_counts_of_taxa_in_align.py

Created by Brant Faircloth on 16 August 2011.
Copyright 2011 Brant C. Faircloth. All rights reserved.

The program iterates through a folder of nexus files and returns the count
of alignments having more than "--percent" of "--taxa" taxa.
"""


import pdb
import os
import sys
import glob
import numpy
import shutil
import argparse
from Bio import AlignIO
from phyluce.helpers import is_dir



def get_args():
    parser = argparse.ArgumentParser(description='Match UCE probes to assembled contigs and store the data')
    parser.add_argument('nexus', help='The directory containing the nexus files', type=is_dir)
    parser.add_argument('taxa', type=int, help='The number of taxa expected')
    parser.add_argument('output', type=is_dir, help='The output dir in which to store copies of the alignments')
    parser.add_argument('--percent', dest = 'percent', type=float, default = 0.5, help='The percent of taxa to require')
    return parser.parse_args()


def get_files(input_dir):
    return glob.glob(os.path.join(os.path.expanduser(input_dir), '*.nex'))

def main():
    args = get_args()
    # iterate through all the files to determine the longest alignment
    files = get_files(args.nexus)
    #align_lengths = [[AlignIO.read(f, 'nexus').get_alignment_length(),os.path.split(f)[1]] \
    #                        for f in files]
    counts = []
    frac = []
    for f in files:
        aln = AlignIO.read(f, 'nexus')
        proportion = round(args.percent * args.taxa)
        if len(aln) < int(proportion):
            pass
        else:
            counts.append(1)
            frac.append(len(aln)/float(args.taxa))
            shutil.copyfile(f,os.path.join(args.output, os.path.basename(f)))
    print "Copied {0} alignments containing an proportion of {1} taxa".format(sum(counts), round(sum(frac)/len(frac), 2))
    #pdb.set_trace()

if __name__ == '__main__':
    main()
