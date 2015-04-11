#!/usr/bin/env python
# encoding: utf-8
"""
File: extract_taxon_data_from_fasta_file.py
Author: Brant Faircloth

Created by Brant Faircloth on 14 August 2012 11:08 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: 

"""


import os
import glob
import argparse
from Bio import SeqIO
from phyluce.helpers import FullPaths
#import pdb


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Extract sequence of a given taxa from alignments""")
    parser.add_argument(
            "fasta",
            action=FullPaths,
            help="""The file/directory of alignments"""
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
    return parser.parse_args()


def parse_fasta_file(f, taxon, outfile):
    for seq in SeqIO.parse(f, 'fasta'):
        if taxon in seq.id:
            outfile.write(seq.format('fasta'))


def main():
    args = get_args()
    if os.path.isdir(args.fasta):
        for count, f in enumerate(glob.glob(os.path.join(args.fasta, "*.fa*"))):
            print count
            parse_fasta_file(f, args.taxon, args.output)
    elif os.path.isfile(args.fasta):
        parse_fasta_file(args.fasta, args.taxon, args.output)
    args.output.close()


if __name__ == '__main__':
    main()
