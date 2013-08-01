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
from Bio.Nexus import Nexus
from phyluce.helpers import is_dir, FullPaths

#import pdb

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
    parser.add_argument(
        "--nexus",
        action="store_true",
        default=False,
        help="""Export as NEXUS format""",
    )
    return parser.parse_args()


def main():
    args = get_args()
    print "Reading files..."
    nexus_files = glob.glob(os.path.join(args.input, '*.nex*'))
    data = [(fname, Nexus.Nexus(fname)) for fname in nexus_files]
    print "Concatenating files..."
    concatenated = Nexus.combine(data)
    if not args.nexus:
        print "Writing to PHYLIP format..."
        concatenated.export_phylip(args.output)
    else:
        print "Writing to NEXUS format..."
        concatenated.write_nexus_data(args.output)


if __name__ == '__main__':
    main()
