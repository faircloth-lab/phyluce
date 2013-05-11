#!/usr/bin/env python
# encoding: utf-8
"""
File: move_align_by_conf_file.py
Author: Brant Faircloth

Created by Brant Faircloth on 08 June 2012 11:06 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: 

"""

import os
import glob
import shutil
import argparse
import ConfigParser
from phyluce.helpers import is_dir, FullPaths

import pdb


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Move files if they are in (or with --opposite, not in) a config file""")
    parser.add_argument(
            "conf",
            help="""The configuration file giving locus names"""
        )
    parser.add_argument(
            "input",
            action=FullPaths,
            type=is_dir,
            help="""The input directory for the alignments"""
        )
    parser.add_argument(
            "output",
            action=FullPaths,
            type=is_dir,
            help="""The output directory for the alignments"""
        )
    parser.add_argument(
            "--sections",
            nargs='+',
            default=None,
            help="""The sections of files to use as a filter (defaults to all)""",
        )
    parser.add_argument(
            "--opposite",
            action="store_true",
            default=False,
            help="""Move alignments based on what is NOT in the conf file""",
        )
    parser.add_argument(
            "--extension",
            type=str,
            default="nex",
            help="""The extension of the files to move""",
        )
    return parser.parse_args()


def main():
    args = get_args()
    conf = ConfigParser.ConfigParser()
    conf.read(args.conf)
    if not args.sections:
        args.sections = conf.sections()
    #pdb.set_trace()
    temp_items = []
    for section in args.sections:
        temp_items.extend([i[0] for i in conf.items(section)])
    items = set(temp_items)
    for src in glob.glob(os.path.join(args.input, '*.{}*'.format(args.extension))):
        basename = os.path.basename(src)
        if not args.opposite:
            if os.path.splitext(basename)[0] in items:
                shutil.copyfile(src, os.path.join(args.output, basename))
            else:
                pass
        else:
            if os.path.splitext(basename)[0] not in items:
                shutil.copyfile(src, os.path.join(args.output, basename))
            else:
                pass

if __name__ == '__main__':
    main()