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
from phyluce.helpers import is_dir, FullPaths

import pdb


def get_args():
    parser = argparse.ArgumentParser(description="""Match UCE probes
                to assembled contigs and store the data""")
    parser.add_argument('nexus',
        type=is_dir,
        help="""The directory containing the nexus files""")
    parser.add_argument('--output',
        type=is_dir,
        action=FullPaths,
        help="""The directory to copy good alignments to""")
    return parser.parse_args()


def get_files(input_dir):
    return glob.glob(os.path.join(os.path.expanduser(input_dir), '*.nex'))


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
    files = get_files(args.nexus)
    # compile our regexes once
    n_bases = re.compile("N{3,}")
    for f in files:
        aln = AlignIO.read(f, 'nexus')
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
