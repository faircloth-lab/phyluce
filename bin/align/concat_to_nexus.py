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
            type=is_dir,
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
    aln = AlignIO.read(args.nexus, 'nexus')
    print "Writing output"
    for k in sorted(partitions.iterkeys()):
        #sys.stdout.write('.')
        #sys.stdout.flush()
        try:
            start, end = partitions[k][0], partitions[k][-1]
            #pdb.set_trace()
            #print start, end
            temp = aln[:, start:end]
            outf = open(os.path.join(args.output, "{}.nex".format(k)), 'w')
            outf.write(temp.format('nexus'))
        except IndexError:
            print "Died on partition {0}".format(k)
        outf.close()

# charset AC010974= 1320460-1320459;
# charset AC069154= 1324339-1324338;
if __name__ == '__main__':
    main()
