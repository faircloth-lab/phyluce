#!/usr/bin/env python
# encoding: utf-8
"""
File: get_genome_sequences_from_bed.py
Author: Brant Faircloth

Created by Brant Faircloth on 09 June 2013 11:06 PDT (-0700)
Copyright (c) 2013 Brant C. Faircloth. All rights reserved.

Description:

"""

import argparse

from bx.seq import twobit
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phyluce.helpers import FullPaths

import pdb


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
        description="""Given an input BED file for one taxon, slice out genomic sequence given coordinates"""
    )
    parser.add_argument(
        "--bed",
        required=True,
        help="""The input BED file"""
    )
    parser.add_argument(
        "--twobit",
        required=True,
        action=FullPaths,
        help="""The input genome to slice in UCSC 2bit format"""
    )
    parser.add_argument(
        "--output",
        required=True,
        help="""The output fasta file"""
    )
    parser.add_argument(
        "--filter-mask",
        dest='mask',
        type=float,
        default=0.25,
        help="""Filter strings with > X frequency of masked bases""",
    )
    parser.add_argument(
        "--max-n",
        type=int,
        default=0,
        help="""The maximum number of ambiguous bases ('N') to accept""",
    )
    parser.add_argument(
        "--buffer-to",
        type=int,
        default=None,
        help="""The length to which to buffer the extracted sequences""",
    )
    return parser.parse_args()


def create_sequence_object(cnt, sequence, chromo, start, end):
    name = "slice_{} |{}:{}-{}".format(cnt, chromo, start, end)
    return SeqRecord(Seq(sequence), id=name, name='', description='')


def sequence_is_masked(mask, sequence, cnt):
    if mask == None:
        return False
    else:
        masked_bases = float(sum([1 for i in sequence if i.islower()]))
        if (masked_bases/len(sequence)) > mask:
            return True
        else:
            return False

def n_count(sequence):
    return sequence.upper().count('N')


def main():
    args = get_args()
    tb = twobit.TwoBitFile(file(args.twobit))
    filtered = 0
    kept = 0
    with open(args.output, 'w') as outf:
        for cnt, line in enumerate(open(args.bed, 'rU')):
            ls = [int(i) if i.isdigit() else i for i in line.strip().split('\t')]
            chromo, start, end = ls
            if args.buffer_to is not None:
                length = abs(end - start)
                delta = args.buffer_to - length
                if delta > 0:
                    if delta % 2 != 0:
                        delta += 1
                    start -= delta / 2
                    end += delta / 2
            sequence = tb[str(chromo)][start:end]
            if len(sequence) >= args.buffer_to and n_count(sequence) <= args.max_n and not sequence_is_masked(args.mask, sequence, cnt):
                seq = create_sequence_object(cnt, sequence, chromo, start, end)
                outf.write(seq.format('fasta'))
                kept += 1
            else:
                filtered += 1
    print "Screened {} sequences.  Filtered {} < 160 bp or with > {}% masked bases or > {} N-bases. Kept {}.".format(
        cnt + 1,
        filtered,
        args.mask * 100,
        args.max_n,
        kept
    )


if __name__ == '__main__':
    main()
