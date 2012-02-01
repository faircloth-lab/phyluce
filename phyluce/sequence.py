#!/usr/bin/env python
# encoding: utf-8

"""
sequence.py

Created by Brant Faircloth on 04 May 2010 21:34 PDT (-0700).
Copyright (c) 2010 Brant C. Faircloth. All rights reserved.
"""

import pdb
import os
import sys
import string
import hashlib
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord


def reverse_complement(seq):
    '''Return reverse complement of seq'''
    bases = string.maketrans('AGCTagct','TCGAtcga')
    # translate it, reverse, return
    return seq.translate(bases)[::-1]

def formatter(options, tb, handle, sp, name1, strand1, zstart1, end1, name2, strand2, zstart2, end2, debug=False):
    if debug:
        pdb.set_trace()
    preceding   = zstart1 - zstart2 - options.length
    # make sure we cannot slice sequence that does not exist
    if preceding < 0:
        preceding = 0
    if options.uce:
        following = end1 + options.length
    else:
        following   = end1 + (120 - end2) + options.length
    entire      = tb[name1][preceding:following]
    if strand1 == '+' and strand2 == '+':
        dna     = entire
    elif strand1 == '+' and strand2 == '-':
        dna     = reverse_complement(entire)
    # this is the position of the sequence were drawing
    thisSeqString = '%s:%s-%s|%s:%s-%s' % (name1, preceding, following, name1, zstart1, end1)
    longName = '%s|match to bases %s to %s of|%s' % (thisSeqString, zstart2, end2, name2)
    seq = Seq(dna, IUPAC.unambiguous_dna)
    record = SeqRecord(seq)
    # get unique record id as hash of UCSC map position
    record.id = '%s_%s' % (sp, hashlib.md5(thisSeqString).hexdigest()[-6:])
    record.name = record.id
    record.description = longName
    # go ahead and write it
    handle.write(record.format('fasta'))


if __name__ == '__main__':
    main()