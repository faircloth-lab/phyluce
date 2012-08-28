#!/usr/bin/env python
# encoding: utf-8
"""
File: extract_taxon_data_from_alignments.py
Author: Brant Faircloth

Created by Brant Faircloth on 02 June 2012 14:06 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: Given a directory of alignments, extract those
sequence data for a single taxon and output in fasta format.

"""

import os
import glob
import argparse
from Bio import AlignIO
from phyluce.helpers import is_dir, FullPaths, get_file_extensions

import pdb


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Extract sequence of a given taxa from alignments""")
    parser.add_argument(
            "alignments",
            type=is_dir,
            action=FullPaths,
            help="""The directory of alignments"""
        )
    parser.add_argument(
            "taxon",
            type=str,
            help="""The taxon to extract"""
        )
    parser.add_argument(
            "output",
            type=argparse.FileType('w'),
            help="""The output file"""
        )
    parser.add_argument(
            "--input-format",
            dest='input_format',
            choices=['nexus', 'newick', 'fasta', 'phylip'],
            default='nexus',
            help="""The input format of the alignments""",
        )
    return parser.parse_args()


def main():
    args = get_args()
    alignments = []
    for ftype in get_file_extensions(args.input_format):
        alignments.extend(glob.glob(os.path.join(args.alignments, "*{}".format(ftype))))
    for count, f in enumerate(alignments):
        aln = AlignIO.read(f, args.input_format)
        for taxon in aln:
            if taxon.id == args.taxon:
                seq = str(taxon.seq).replace('-', '')
                locus = os.path.splitext(os.path.basename(f))[0]
                if not len(seq) == 0:
                    args.output.write(">{0}\n{1}\n".format(locus, seq))
                else:
                    print locus
    args.output.close()


if __name__ == '__main__':
    main()

