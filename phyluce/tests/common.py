#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2018 Brant Faircloth || http://faircloth-lab.org/

All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 06 July 2018 15:27 CDT (-0500)
"""

from Bio import SeqIO


def get_contig_lengths_and_counts(contigs):
    with open(contigs) as contig_file:
        contig_count = 0
        contig_length = 0
        for seq in SeqIO.parse(contig_file, 'fasta'):
            contig_count += 1
            contig_length += len(seq)
    return contig_count, contig_length
