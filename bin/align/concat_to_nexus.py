#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2015 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 09 April 2015 14:17 CDT (-0500)
"""

import os
import argparse
from Bio import AlignIO
from Bio.Nexus import Nexus
from phyluce.helpers import FullPaths, is_file, CreateDir

#import pdb


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Split a concatenated NEXUS file into component loci using the charsets values""")
    parser.add_argument(
            "--nexus",
            required=True,
            type=is_file,
            action=FullPaths,
            help="""The input NEXUS file"""
        )
    parser.add_argument(
            "--output",
            action=CreateDir,
            required=True,
            help="""The output directory in which to store alignments""",
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


if __name__ == '__main__':
    main()
