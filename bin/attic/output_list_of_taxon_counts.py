#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2014 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 19 December 2014 11:13 CST (-0600)
"""


import os
import glob
import argparse
import multiprocessing
from Bio import SeqIO

from phyluce.helpers import is_dir, FullPaths

import pdb


def get_args():
    parser = argparse.ArgumentParser(
        description="""Compute summary statistics for alignments in parallel""",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--fastas",
        required=True,
        type=is_dir,
        action=FullPaths,
        help="""The directory containing fastas to checked."""
    )
    parser.add_argument(
        "--input-format",
        dest="input_format",
        choices=["fasta"],
        default="fasta",
        help="""The input file format.""",
    )
    parser.add_argument(
        "--min-count",
        type=int,
        default=3,
        help="""The min count of taxa allowed in a fasta file"""
    )
    parser.add_argument(
        "--cores",
        type=int,
        default=1,
        help="""Process alignments in parallel using --cores for alignment. """ +
        """This is the number of PHYSICAL CPUs."""
    )
    return parser.parse_args()


def get_fasta_files(input_dir, input_format):
    fastas = []
    for ftype in ('.fasta', '.fsa', '.fa'):
        fastas.extend(glob.glob(os.path.join(input_dir, "*{}".format(ftype))))
    return fastas


def count_taxa(work):
    file, format = work
    with open(file, "rU") as infile:
        counts = [1 for seq in SeqIO.parse(infile, format)]
    return file, sum(counts)

def main():
    args = get_args()
    files = get_fasta_files(args.fastas, args.input_format)
    work = [[file, args.input_format] for file in files]
    if args.cores > 1:
        assert args.cores <= multiprocessing.cpu_count(), "You've specified more cores than you have"
        pool = multiprocessing.Pool(args.cores)
        summary = pool.map(count_taxa, work)
    else:
        summary = map(count_taxa, work)
    for file, count in summary:
        if count < args.min_count:
            print "{}".format(os.path.basename(file))

if __name__ == '__main__':
    main()
