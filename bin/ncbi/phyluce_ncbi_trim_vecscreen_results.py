#!/usr/bin/env python
# encoding: utf-8
"""
File: trim_vecscreen_results.py
Author: Brant Faircloth

Created by Brant Faircloth on 18 December 2012 16:12 PST (-0800)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: Trim vector contamination from data using a blast
output file.  Run blast like so:

blastn -task blastn -db UniVec_core -query test.fsa \
    -evalue 1 -gapopen 3 -gapextend 3 -word_size 11 \
    -reward 1 -penalty -5 -out blast.out -num_threads 4 \
    -dust yes -searchsp 1750000000000 -soft_masking true \
    -outfmt 6

Run this program like so:

python ~/phyluce/bin/ncbi/trim_vecscreen_results.py \
    blast.out \
    output.fsa \
    output-univec-trimmed.fsa

"""

import argparse
from seqtools.sequence import FastaReader, FastaWriter
from collections import namedtuple, defaultdict

#import pdb


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Program description""")
    parser.add_argument(
            "blast",
            help="""The input blast results (-outfmt 6)"""
        )
    parser.add_argument(
            "fasta",
            help="""The fasta file to trim"""
        )
    parser.add_argument(
            "output",
            help="""The output file to save"""
        )
    return parser.parse_args()


def main():
    args = get_args()
    hits = defaultdict(list)
    for line in open(args.blast, 'rU'):
        BlastRecord = namedtuple(
            "BlastRecord",
            "query_id,subject_id,perc_identity,aln_length,mismatch_count,gap_open_count,query_start,query_end,subject_start,subject_end,e_val,bit_score"
        )
        rec = BlastRecord._make(line.strip().split("\t"))
        hits[rec.query_id].append(rec)
    outf = FastaWriter(args.output)
    for sequence in FastaReader(args.fasta):
        name = sequence.identifier.split(' ')[0].strip('>')
        if name not in hits.keys():
            outf.write(sequence)
        else:
            oldlen = len(sequence)
            front = []
            back = []
            for hit in hits[name]:
                if int(hit.query_start) < len(sequence) / 2:
                    front.append(int(hit.query_end))
                else:
                    back.append(int(hit.query_start))
            if front:
                # get max start pos of matches at start of sequence
                fmx = max(front)
            else:
                fmx = 0
            if back:
                # get min start pos of matches at end of sequence
                bmn = min(back) - 1
            else:
                bmn = len(sequence)
            # slice sequence
            sequence = sequence.slice(fmx, bmn)
            outf.write(sequence)
            print "Trimmed sequence {0}, start={1}, end={2}".format(name, fmx, bmn - oldlen)


if __name__ == '__main__':
    main()