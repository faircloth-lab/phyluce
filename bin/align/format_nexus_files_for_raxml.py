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

from phyluce.helpers import FullPaths, CreateDir, is_dir
from phyluce.log import setup_logging

#import pdb

def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
        description="""Convert individual nexus files to concatenated phylip format""",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--alignments",
        required=True,
        type=is_dir,
        action=FullPaths,
        help="""The directory containing alignments to concatenate (NEXUS-ONLY)."""
    )
    parser.add_argument(
        "--output",
        required=True,
        action=CreateDir,
        help="""The output file for the concatenated phylip data""",
    )
    parser.add_argument(
        "--nexus",
        action="store_true",
        default=False,
        help="""Export as NEXUS format""",
    )
    parser.add_argument(
        "--charsets",
        action="store_true",
        default=False,
        help="""Add charsets to phylip file""",
    )
    parser.add_argument(
        "--verbosity",
        type=str,
        choices=["INFO", "WARN", "CRITICAL"],
        default="INFO",
        help="""The logging level to use."""
    )
    parser.add_argument(
        "--log-path",
        action=FullPaths,
        type=is_dir,
        default=None,
        help="""The path to a directory to hold logs."""
    )
    return parser.parse_args()


def main():
    args = get_args()
    # setup logging
    log, my_name = setup_logging(args)
    # read alignments
    log.info("Reading input alignments in NEXUS format")
    nexus_files = glob.glob(os.path.join(args.alignments, '*.nex*'))
    data = [(os.path.basename(fname), Nexus.Nexus(fname)) for fname in nexus_files]
    log.info("Concatenating files")
    concatenated = Nexus.combine(data)
    if not args.nexus:
        concat_file = os.path.join(args.output, os.path.basename(args.alignments) + ".phylip")
        if args.charsets:
            sets = concatenated.append_sets()
            charset_file = os.path.join(args.output, os.path.basename(args.alignments) + ".charsets")
            log.info("Writing charsets to {}".format(
                charset_file
            ))
            with open(charset_file, 'w') as outf:
                outf.write(sets)
        log.info("Writing concatenated PHYLIP alignment to {}".format(concat_file))
        concatenated.export_phylip(concat_file)
    else:
        concat_file = os.path.join(args.output, os.path.basename(args.alignments) + ".nexus")
        if args.charsets:
            log.info("Writing concatenated alignment to NEXUS format (with charsets)")
            concatenated.write_nexus_data(concat_file)
        else:
            log.info("Writing concatenated alignment to NEXUS format (without charsets)")
            concatenated.write_nexus_data(concat_file, append_sets=False)
    # end
    text = " Completed {} ".format(my_name)
    log.info(text.center(65, "="))


if __name__ == '__main__':
    main()
