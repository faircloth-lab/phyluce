#!/usr/bin/env python
# encoding: utf-8
"""
File: get_fish_probe_overlaps.py
Author: Brant Faircloth

Created by Brant Faircloth on 04 October 2012 13:10 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: 

"""

import argparse
from Bio import SeqIO
from phyluce.lastz import Reader
from collections import defaultdict

import pdb

def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Program description""")
    parser.add_argument(
            "lastz",
            help="""A lastz file mapping first to second set of fish probes"""
        )
    parser.add_argument(
            "fasta1",
            help="""The fasta file of one set of probes"""
        )
    parser.add_argument(
            "fasta2",
            help="""The fasta file of the 2nd set of probes"""
        )
    parser.add_argument(
            "output",
            help="""The output fasta file"""
        )
    return parser.parse_args()

def main():
    args = get_args()
    lastz = Reader(args.lastz)
    identical = defaultdict(list)
    for row in lastz:
        n1 = row.name1.split('|')[0]
        n2 = row.name2.split('|')[0]
        if n1 == n2:
            span1 = row.name1.split('|')[-1]
            span2 = row.name2.split('|')[-1]
            if span1 == span2:
                identical[n1].append(span1)
    new_fasta = open(args.output, 'w')
    for record in SeqIO.parse(open(args.fasta1, 'rU'), 'fasta'):
        rs = record.id.split('|')
        name, positions = rs[0], rs[-1]
        #pdb.set_trace()
        if positions in identical[name]:
            record.description = '|'.join([rs[0], rs[-1], "both"])
        else:
            record.description = '|'.join([rs[0], rs[-1], "500"])
        new_fasta.write(">{0}\n{1}\n".format(record.description, str(record.seq)))
    for record in SeqIO.parse(open(args.fasta2, 'rU'), 'fasta'):
        rs = record.id.split('|')
        name, positions = rs[0], rs[-1]
        if positions in identical[name]:
            pass
        else:
            record.description = '|'.join([rs[0], rs[-1], "1000"])
            new_fasta.write(">{0}\n{1}\n".format(record.description, str(record.seq)))
    new_fasta.close()


if __name__ == '__main__':
    main()