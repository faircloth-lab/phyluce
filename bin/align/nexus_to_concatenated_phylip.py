#!/usr/bin/env python
# encoding: utf-8
"""
File: nexus_to_concatenated_phylip.py
Author: Brant Faircloth

Created by Brant Faircloth on 21 April 2012 16:04 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: Convert individual files to a large concantenated phylip file

"""

import os
import glob
import argparse
import tempfile

from Bio import AlignIO
from Bio.Nexus import Nexus
from phyluce.helpers import is_dir, FullPaths


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Convert individual nexus files to concatenated phylip format""")
    parser.add_argument(
            "input",
            type=is_dir,
            action=FullPaths,
            help="""The input directory of nexus files"""
        )
    parser.add_argument(
            "output",
            action=FullPaths,
            help="""The output file for the concatenated phylip data""",
        )
    return parser.parse_args()


def main():
    args = get_args()
    print "Reading files..."
    nexus_files = glob.glob(os.path.join(args.input, '*.nex*'))
    data = [(fname, Nexus.Nexus(fname)) for fname in nexus_files]
    print "Concatenating files..."
    concatenated = Nexus.combine(data)
    #print "Writing temp nexus..."
    #fd, temp = tempfile.mkstemp(suffix='.nexus')
    #concatenated.write_nexus_data(filename=os.fdopen(fd, 'w'))
    print "Writing to phylip..."
    concatenated.export_phylip(args.output)
    #aln = AlignIO.parse(temp, "nexus")
    #AlignIO.write(aln, args.output, "phylip")
    #print "Cleaning up..."
    #os.remove(temp)


if __name__ == '__main__':
    main()
