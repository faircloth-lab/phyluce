#!/usr/bin/env python
# encoding: utf-8
"""
File: get_informative_sites.py
Author: Brant Faircloth

Created by Brant Faircloth on 07 August 2012 21:08 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: 

"""

import os
import glob
import argparse
import multiprocessing
from Bio import AlignIO
from collections import Counter
from phyluce.helpers import is_dir, FullPaths, get_file_extensions


import pdb


def get_args():
    parser = argparse.ArgumentParser(
            description="""Count the number of informative sites in a given set of alignment"""
        )
    parser.add_argument(
            'input',
            type=is_dir,
            action=FullPaths,
            help="""The directory containing the alignment files"""
        )
    parser.add_argument(
            '--output',
            type=str,
            default=None,
            help="""The output filename"""
        )
    parser.add_argument(
            "--input-format",
            dest="input_format",
            choices=['fasta', 'nexus', 'phylip', 'clustal', 'emboss', 'stockholm'],
            default='fasta',
            help="""The input alignment format""",
        )
    parser.add_argument(
            "--cores",
            type=int,
            default=1,
            help="""The number of cores to use.""",
        )
    return parser.parse_args()


def get_files(input_dir, input_format):
    alignments = []
    for ftype in get_file_extensions(input_format):
        alignments.extend(glob.glob(os.path.join(input_dir, "*{}".format(ftype))))
    return alignments


def get_informative_sites(count):
    # remove gaps
    del count['-']
    # remove N
    del count['N']
    # remove ?
    del count['?']
    sufficient_sites = len(count)
    if sufficient_sites >= 2:
        sufficient_sequences = sum([1 for i in count.values() if i >= 2])
        if sufficient_sequences >= 2:
            return True
    return False

def get_differences(count):
    # remove gaps
    del count['-']
    # remove N
    del count['N']
    # remove ?
    del count['?']
    sufficient_sites = len(count)
    if sufficient_sites >= 2:
        return True
    else:
        return False

def worker(work):
    args, f = work
    aln = AlignIO.read(f, args.input_format)
    name = os.path.basename(f)
    informative_sites = []
    differences = []
    for idx in xrange(aln.get_alignment_length()):
        col = aln[:, idx]
        count = Counter(col)
        if get_informative_sites(count):
            informative_sites.append(1)
        else:
            informative_sites.append(0)
        if get_differences(count):
            differences.append(1)
        else:
            differences.append(0)
    return (name, aln.get_alignment_length(), sum(informative_sites), sum(differences))


def main():
    args = get_args()
    work = [(args, f) for f in get_files(args.input, args.input_format)]
    if args.cores <= 1:
        results = map(worker, work)
    elif args.cores > 1:
        pool = multiprocessing.Pool(args.cores)
        results = pool.map(worker, work)
    if args.output:
        outf = open(args.output, 'w')
        outf.write("locus,length,informative_sites,differences\n")
    else:
        print "locus\tlength\tinformative_sites\tdifferences"
    for locus in results:
        if not args.output:
            print "{0}\t{1}\t{2}\t{3}".format(locus[0], locus[1], locus[2], locus[3])
        else:
            outf.write("{0},{1},{2},{3}\n".format(locus[0], locus[1], locus[2], locus[3]))

if __name__ == '__main__':
    main()
