"""
File: interleave_reads.py
Author: Brant Faircloth

Created by Brant Faircloth on 03 March 2012 14:03 PST (-0800)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description:  Given two fastq files of reads, merge those two files
into a single "shuffled" or "interleaved" file of reads.

"""

import sys
import argparse
from itertools import izip
from seqtools.sequence import fastq

def get_args():
    parser = argparse.ArgumentParser(description="Interleave two split, " + \
        "paired-end fastq file into one file")
    parser.add_argument("read1", help="The output read1 FASTQ file name")
    parser.add_argument("read2", help="The output read2 FASTQ file name")
    parser.add_argument("output", help="The output fastq file")
    args = parser.parse_args()
    return args


def main():
    args = get_args()
    outfile = fastq.FasterFastqWriter(args.output)
    read1 = fastq.FasterFastqReader(args.read1)
    read2 = fastq.FasterFastqReader(args.read2)
    rc = 0
    sys.stdout.write("Interleaving reads (1 dot = 10,000 pairs): ")
    for r1, r2 in izip(read1, read2):
        assert r1[0].split(" ")[0] == r2[0].split(" ")[0], \
                "Read FASTQ headers mismatch."
        outfile.write(r1)
        outfile.write(r2)
        if rc != 0 and rc % 10000 == 0:
            sys.stdout.write(".")
            sys.stdout.flush()
        rc += 1
    outfile.close()
    read1.close()
    read2.close()

if __name__ == '__main__':
    main()
