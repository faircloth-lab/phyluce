#!/usr/bin/env python
# encoding: utf-8

"""
get_align_summary_data.py

Created by Brant Faircloth on 16 August 2011.
Copyright 2011 Brant C. Faircloth. All rights reserved.

The program iterates through a folder of nexus files and returns the count
of alignments having more than "--percent" of "--taxa" taxa.
"""


import os
import sys
import glob
import numpy
import shutil
import argparse
from Bio import AlignIO
from collections import Counter, defaultdict
from phyluce.helpers import is_dir, FullPaths, get_file_extensions

#import pdb


def get_args():
    parser = argparse.ArgumentParser(
            description="""Compute summary parameters for alignments"""
        )
    parser.add_argument(
            'input',
            type=is_dir,
            action=FullPaths,
            help='The directory containing the alignment files'
        )
    parser.add_argument(
            '--min-taxa',
            type=int,
            help='The minimum number of taxa to count'
        )
    parser.add_argument(
            "--input-format",
            dest="input_format",
            choices=['fasta', 'nexus', 'phylip', 'clustal', 'emboss', 'stockholm'],
            default='fasta',
            help="""The input alignment format""",
        )
    return parser.parse_args()


def get_files(input_dir, input_format):
    alignments = []
    for ftype in get_file_extensions(input_format):
        alignments.extend(glob.glob(os.path.join(input_dir, "*{}".format(ftype))))
    return alignments


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
    files = get_files(args.input, args.input_format)
    counts = []
    lengths = []
    trimmed = []
    bases = defaultdict(list)
    for f in files:
        try:
            aln = AlignIO.read(f, args.input_format)
            l = aln.get_alignment_length()
            lengths.append(l)
            if l < 100:
                print "{0} is < 100 bp long".format(os.path.basename(f))
            for col in xrange(len(aln[0])):
                for base in aln[:, col]:
                    bases[base.upper()].append(1)
            if args.min_taxa:
                if len(aln) >= args.min_taxa:
                    counts.append(len(aln))
                else:
                    print "{0} has fewer than {1} taxa".format(f, args.min_taxa)
            else:
                counts.append(len(aln))
            for read in aln:
                #pdb.set_trace()
                left = len(read.seq) - len(str(read.seq).lstrip('-'))
                right = len(read.seq) - len(str(read.seq).rstrip('-'))
                trimmed.append(left + right)
        except ValueError, e:
            if e.message == 'No records found in handle':
                print 'No records found in {0}'.format(os.path.basename(f))
            else:
                raise ValueError('Something is wrong with alignment {0}'.format(os.path.basename(f)))
    compute_lengths(lengths)
    compute_taxa(counts)
    compute_bases(bases, trimmed)


if __name__ == '__main__':
    main()
