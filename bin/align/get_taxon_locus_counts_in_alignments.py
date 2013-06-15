#!/usr/bin/env python
# encoding: utf-8
"""
File: get_taxon_locus_counts_in_alignments.py
Author: Brant Faircloth

Created by Brant Faircloth on 15 June 2013 10:06 PDT (-0700)
Copyright (c) 2013 Brant C. Faircloth. All rights reserved.

Description:

"""

import os
import re
import glob
import argparse
from Bio import AlignIO
from phyluce.helpers import is_dir, FullPaths, get_file_extensions
from collections import Counter

import pdb


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
        description="""Given an input directory of alignments, determine the number of alignments per taxon"""
    )
    parser.add_argument(
        "--input",
        required=True,
        action=FullPaths,
        type=is_dir,
        help="""The input directory of alignment files"""
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


def main():
    args = get_args()
    # iterate through all the files to determine the longest alignment
    files = get_files(args.input, args.input_format)
    count = Counter()
    for f in files:
        aln = AlignIO.read(f, args.input_format)
        for seq in aln:
            if len(set(str(seq.seq))) > 1:
                count[seq.id] += 1
    print "taxon,count"
    for taxon, cnt in count.iteritems():
        print "{},{}".format(taxon, cnt)


if __name__ == '__main__':
    main()
