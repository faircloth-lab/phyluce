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

from collections import Counter, defaultdict


def get_args():
    parser = argparse.ArgumentParser(description='Match UCE probes to assembled contigs and store the data')
    parser.add_argument('nexus', help='The directory containing the nexus files', type=is_dir)
    parser.add_argument('--min-taxa', help='''The minimum number of taxa to
            count''', type=int)
    parser.add_argument('--count-bases', dest='count_bases', help='''The minimum number of taxa to
            count''', action='store_true', default=False)
    return parser.parse_args()

def get_files(input_dir):
    return glob.glob(os.path.join(os.path.expanduser(input_dir), '*.nex'))

def main():
    args = get_args()
    # iterate through all the files to determine the longest alignment
    files = get_files(args.nexus)
    counts = []
    bases = defaultdict(list)
    for f in files:
        aln = AlignIO.read(f, 'nexus')
        if args.count_bases:
            for col in xrange(len(aln[0])):
                for base in aln[:,col]:
                    bases[base].append(1)
        if args.min_taxa:
            if len(aln) >= args.min_taxa:
                counts.append(len(aln))
        else:
             counts.append(len(aln))
    print "Average: {}".format(sum(counts)/float(len(counts)))
    ci = 1.96 * numpy.std(numpy.array(counts), ddof=1)/numpy.sqrt(len(counts))
    base_sum = {base:sum(bases[base]) for base in bases}
    print "95 CI: {}".format(ci)
    print "(min): {}".format(min(counts))
    print "(max): {}".format(max(counts))
    print "Count {{Num taxa:alns with Num taxa}}: {}".format(dict(Counter(counts)))
    print "Bases: {}".format(base_sum)

if __name__ == '__main__':
    main()
