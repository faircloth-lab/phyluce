#!/usr/bin/env python
# encoding: utf-8
"""
File: screen_alignments_for_problems.py
Author: Brant Faircloth

Created by Brant Faircloth on 07 March 2012 15:03 PST (-0800)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: 

"""

import os
import re
import glob
import shutil
import argparse
from Bio import AlignIO
from phyluce.helpers import is_dir, FullPaths, get_file_extensions

import pdb


def get_args():
    parser = argparse.ArgumentParser(description="""Screen a given set of alignments
        for problem bases""")
    parser.add_argument('input',
        type=is_dir,
        action=FullPaths,
        help="""The directory containing the nexus files""")
    parser.add_argument(
            "--input-format",
            dest="input_format",
            choices=['fasta', 'nexus', 'phylip', 'clustal', 'emboss', 'stockholm'],
            default='fasta',
            help="""The input alignment format""",
        )
    parser.add_argument('--output',
        type=is_dir,
        action=FullPaths,
        help="""The directory to copy good alignments to""")
    return parser.parse_args()


def get_files(input_dir, input_format):
    alignments = []
    for ftype in get_file_extensions(input_format):
        alignments.extend(glob.glob(os.path.join(input_dir, "*{}".format(ftype))))
    return alignments


def find_multiple_N_bases(regex, aln):
    n = {}
    for seq in list(aln):
        result = regex.findall(str(seq.seq))
        if result:
            ncount = [len(i) for i in result]
            n[seq.name] = ncount
    return n


def find_any_X_bases(f, aln):
    x_bases = False
    for seq in list(aln):
        if 'X' in str(seq.seq) or 'x' in str(seq.seq):
            x_bases = True
            print "X-bases in ", os.path.basename(f)
    return x_bases


def main():
    args = get_args()
    # iterate through all the files to determine the longest alignment
    files = get_files(args.input, args.input_format)
    # compile our regexes once
    n_bases = re.compile("N{3,}")
    for f in files:
        aln = AlignIO.read(f, args.input_format)
        n = find_multiple_N_bases(n_bases, aln)
        x_bases = find_any_X_bases(f, aln)
        for name, count in n.iteritems():
            print "{0}\n\t{1}\n\t{2}".format(f, name, count)
        if args.output and not x_bases:
            #pdb.set_trace()
            fname = os.path.basename(f)
            outpath = os.path.join(args.output, fname)
            shutil.copy(f, outpath)


if __name__ == '__main__':
    main()
