#!/usr/bin/env python
# encoding: utf-8
"""
File: associate_node_name_with_uce_contig.py
Author: Brant Faircloth

Created by Brant Faircloth on 21 November 2012 15:11 PST (-0800)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: 

"""

import os
import sqlite3
import argparse
import textwrap
from phyluce.helpers import is_file, FullPaths
from Bio import SeqIO

import pdb


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Given a probe.matches.sqlite file and a file of contigs,
            associate UCE contig name with node-name from fasta""")
    parser.add_argument(
            "db",
            type=is_file,
            action=FullPaths,
            help="""The probe.matches.sqlite database"""
        )
    parser.add_argument(
            "fasta",
            type=is_file,
            action=FullPaths,
            help="""The fasta file of contigs"""
        )
    parser.add_argument(
            "output",
            action=FullPaths,
            help="""The output file to hold results"""
        )
    parser.add_argument(
            "--concatenate",
            action="store_true",
            default=False,
            help="""Concatenate sequences rather than making them individual""",
        )
    parser.add_argument(
            "--spacer",
            type=int,
            default=2000,
            help="""The default N-base spacer to insert between concatenated fastas"""
        )
    return parser.parse_args()


def quickwrap(string, length=60):
    return [string[i:i + length] for i in range(0, len(string), length)]


def main():
    args = get_args()
    outf = open(args.output, 'w')
    if not args.concatenate:
        conn = sqlite3.connect(args.db)
        cur = conn.cursor()
        taxon = os.path.splitext(os.path.basename(args.fasta))[0]
        # get all matching contigs for taxon
        query = "SELECT uce, {0} FROM match_map WHERE {0} IS NOT NULL".format(taxon)
        cur.execute(query)
        results = dict([(record[1].split('(')[0], record[0]) for record in cur.fetchall()])
        for seq in SeqIO.parse(open(args.fasta, 'rU'), 'fasta'):
            name = '_'.join(seq.id.split('_')[:2]).lower()
            new_id = results[name]
            seq.id = new_id
            seq.name = ""
            seq.description = ""
            outf.write(seq.format('fasta'))
    else:
        seqs = []
        outf.write('>concatenated-uce-loci\n')
        for seq in SeqIO.parse(open(args.fasta, 'rU'), 'fasta'):
            seqs.append(str(seq.seq))
        spacer = 'N' * args.spacer
        # insert Ns at front and back of list
        seqs.insert(0, spacer)
        seqs.append(spacer)
        # insert other spacer Ns
        big_string = spacer.join(seqs)
        # wrap the sequence to make it pretty
        wrap_list = quickwrap(big_string)
        outf.write('\n'.join(wrap_list))
    outf.close()


if __name__ == '__main__':
    main()