#!/usr/bin/env python
# encoding: utf-8

"""
get_align_summary_data.py

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
    parser = argparse.ArgumentParser(
            description="""Compute summary parameters for alignments"""
        )
    parser.add_argument(
            'nexus',
            type=is_dir,
            help='The directory containing the nexus files'
        )
    parser.add_argument(
            '--min-taxa',
            type=int,
            help='The minimum number of taxa to count'
        )
    return parser.parse_args()


def get_files(input_dir):
    return glob.glob(os.path.join(os.path.expanduser(input_dir), '*.nex'))


def pretty_printer(result):
    print "{:<25}{:<20}".format(result[0], result[1])


def compute_lengths(lengths):
    print "\nLengths\n-----"
    l = numpy.array(lengths)
    sm = ["Total length(aln)", sum(l)]
    avg = ["Average length(aln)", numpy.mean(l)]
    ci = ["95 CI length(aln)", \
            1.96 * (numpy.std(l, ddof=1) / numpy.sqrt(len(l)))
            ]
    for result in [sm, avg, ci]:
        pretty_printer(result)


def compute_taxa(counts):
    print "\nTaxa\n-----"
    avg = ["Average(taxa)", sum(counts) / float(len(counts))]
    ci = [
                "95 CI(taxa)",
                1.96 * numpy.std(numpy.array(counts), ddof=1) / \
                numpy.sqrt(len(counts))
            ]
    mn = ["min(taxa)", min(counts)]
    mx = ["max(taxa)", max(counts)]
    cnt = ["Count(taxa:# alns)", dict(Counter(counts))]
    for result in [avg, ci, mn, mx, cnt]:
        pretty_printer(result)


def compute_bases(bases, trimmed):
    print "\nBase composition\n-----"
    bssm = {base:sum(bases[base]) for base in bases}
    sm = ["Bases", bssm]
    al = ["Sum(all)", sum(bssm.values())]
    nogpsm = sum([bssm[i] for i in ['A', 'C', 'G', 'T']])
    nogp = ["Sum(nucleotide only)", nogpsm]
    trim = ["Missing data from trim (%)", round(sum(trimmed) / float(sum(bssm.values())) * 100, 2)]
    for result in [sm, al, nogp, trim]:
        pretty_printer(result)


def main():
    args = get_args()
    # iterate through all the files to determine the longest alignment
    files = get_files(args.nexus)
    counts = []
    lengths = []
    trimmed = []
    bases = defaultdict(list)
    for f in files:
        aln = AlignIO.read(f, 'nexus')
        lengths.append(aln.get_alignment_length())
        for col in xrange(len(aln[0])):
            for base in aln[:, col]:
                bases[base.upper()].append(1)
        if args.min_taxa:
            if len(aln) >= args.min_taxa:
                counts.append(len(aln))
        else:
            counts.append(len(aln))
        for read in aln:
            #pdb.set_trace()
            left = len(read.seq) - len(str(read.seq).lstrip('-'))
            right = len(read.seq) - len(str(read.seq).rstrip('-'))
            trimmed.append(left + right)
    compute_lengths(lengths)
    compute_taxa(counts)
    compute_bases(bases, trimmed)


if __name__ == '__main__':
    main()
