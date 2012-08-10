#!/usr/bin/env python
# encoding: utf-8
"""
File: check_alignments_for_missing_data.py
Author: Brant Faircloth

Created by Brant Faircloth on 10 August 2012 11:08 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: 

"""

import os
import glob
import argparse
from Bio import AlignIO
from phyluce.helpers import is_dir, FullPaths, get_file_extensions


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
    for f in files:
        try:
            aln = AlignIO.read(f, args.input_format)
            for taxon in aln:
                if set(taxon.seq) == set('-') or set(taxon.seq) == set('?'):
                    print os.path.basename(f)
        except ValueError, e:
            if e.message == 'No records found in handle':
                print 'No records found in {0}'.format(os.path.basename(f))
            else:
                raise ValueError('Something is wrong with alignment {0}'.format(os.path.basename(f)))

if __name__ == '__main__':
    main()
