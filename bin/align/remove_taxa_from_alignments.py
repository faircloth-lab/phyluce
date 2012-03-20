#!/usr/bin/env python
# encoding: utf-8

"""

remove_taxa_from_alignments.py

Created by Brant Faircloth on 21 September 2011 11:44 PDT (-0700).
Copyright (c) 2011 Brant C. Faircloth. All rights reserved.

PURPOSE:  Remove taxa passed on CLI from a folder of nexus alignment
files.

"""

import pdb
import os
import sys
import glob
import argparse
from Bio import AlignIO
from Bio.Alphabet import IUPAC, Gapped
from Bio.Align.Generic import Alignment

def get_taxa_from_list(string):
    """Convert taxa input as string to a list"""
    try:
        times = [str(i) for i in string.split(',')]
    except:
        raise argparse.ArgumentTypeError("Cannot convert time to list of integers")
    return times

def is_dir(dirname):
    if not os.path.isdir:
        msg = "{0} is not a directory".format(dirname)
        raise argparse.ArgumentTypeError(msg)
    else:
        return dirname

def get_args():
    parser = argparse.ArgumentParser(description='Parse fastq files and drop reads containing Ns.')
    parser.add_argument('input', help="The input directory containing nexus files", type=is_dir)
    parser.add_argument('output', help="The directory in which to store the output files", type=is_dir)
    parser.add_argument('--taxa', help="The list of taxa to prune, as comma separated list", type=get_taxa_from_list)
    return parser.parse_args()

def get_files(input_dir):
    return glob.glob(os.path.join(os.path.expanduser(input_dir), '*.nex'))

def main():
    args = get_args()
    nexus_files = get_files(args.input)
    for count, align_file in enumerate(nexus_files):
        new_align = Alignment(Gapped(IUPAC.unambiguous_dna, "-"))
        for align in AlignIO.parse(align_file, "nexus"):
            for taxon in list(align):
                if taxon.name not in args.taxa:
                    new_align.add_sequence(taxon.name, str(taxon.seq))
        outf = os.path.join(args.output, os.path.basename(align_file))
        AlignIO.write(new_align, open(outf, 'w'), 'nexus')
        print count


if __name__ == '__main__':
    main()
