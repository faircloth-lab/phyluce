#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2014 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 30 March 2014 09:28 PDT (-0700)
"""

import argparse

from bx.seq import twobit
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from phyluce.helpers import FullPaths

import pdb


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
        description="""Given an input locus fasta file, buffer loci to desired length."""
    )
    parser.add_argument(
        "--fasta",
        required=True,
        help="""The input FASTA file"""
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
        "--buffer-length",
        type=float,
        default=180,
        help="""The sequence length to buffer to.""",
    )
    parser.add_argument(
        "--filter-mask",
        dest='mask',
        type=float,
        default=None,
        help="""Filter strings with > X frequency of masked bases""",
    )
    parser.add_argument(
        "--max-n",
        type=int,
        default=0,
        help="""The maximum number of ambiguous bases ('N') to accept""",
    )
    return parser.parse_args()


def get_positions_from_coords(s):
    chromo, coords = s.split(':')
    start, end = [int(i) for i in coords.split('-')]
    return chromo, start, end


def sequence_is_masked(mask, sequence):
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


def create_sequence_object(id, sequence, chromo, start, end):
    name = "{} |{}:{}-{}".format(id, chromo, start, end)
    return SeqRecord(Seq(sequence), id=name, name='', description='')


def main():
    args = get_args()
    tb = twobit.TwoBitFile(file(args.twobit))
    filtered = 0
    kept = 0
    skipped = 0
    with open(args.output, 'w') as outf:
        with open(args.fasta, 'rU') as fasta:
            for record in SeqIO.parse(fasta, 'fasta'):
                chromo, start, end = get_positions_from_coords(record.description.split('|')[1])
                delta = args.buffer_length - (end - start)
                if delta < args.buffer_length:
                    split = int(round(delta/2.))
                    new_start = start - split
                    if new_start < 0:
                        new_start = 0
                    new_end = end + split
                    sequence = tb[str(chromo)][new_start:new_end]
                    if n_count(sequence) <= args.max_n and not sequence_is_masked(args.mask, sequence):
                        seq = create_sequence_object(record.id, sequence, chromo, new_start, new_end)
                        outf.write(seq.format('fasta'))
                        kept += 1
                    else:
                        filtered += 1
                else:
                    outf.write(record.format('fasta'))
                    skipped += 1

    print "Total {} sequences.  Expanded {} and filtered {} with > {}% masked bases or > {} masked bases. Kept {}.".format(
        kept + filtered + skipped,
        kept + filtered,
        filtered,
        args.mask * 100,
        args.max_n,
        kept + skipped
    )

if __name__ == '__main__':
    main()
