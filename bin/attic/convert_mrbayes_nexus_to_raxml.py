#!/usr/bin/env python
# encoding: utf-8
"""
File: concat_to_nexus.py
Author: Brant Faircloth

Created by Brant Faircloth on 26 June 2012 22:06 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: 

"""

import os
import sys
import argparse
from Bio import AlignIO
from Bio.Nexus import Nexus
from phyluce.helpers import is_dir, is_file, FullPaths

import pdb


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Program description""")
    parser.add_argument(
            "nexus",
            type=is_file,
            help="""Help text"""
        )
    parser.add_argument(
            "output",
            action=FullPaths,
            help="""Help text""",
        )
    return parser.parse_args()


def main():
    args = get_args()
    # get the partitions from the nexus file
    print "Getting partition information..."
    aln = Nexus.Nexus()
    aln.read(args.nexus)
    partitions = aln.charsets
    print "\tThere are {0} partitions".format(len(partitions))
    print "Parsing alignment..."
    partition_file = os.path.splitext(args.output)[0] + '.partitions'
    outf = open(partition_file, 'w')
    print "Writing partitions..."
    partitions = {}
    for name, bases in aln.charsets.iteritems():
        partitions[bases[0]] = "DNA, {0} = {1}-{2}".format(name, bases[0] + 1, bases[-1] + 1)
    for name in sorted(partitions.keys()):
        outf.write("{0}\n".format(partitions[name]))
    outf.close()
    print "Writing alignment..."
    aln.export_phylip(args.output)

if __name__ == '__main__':
    main()
