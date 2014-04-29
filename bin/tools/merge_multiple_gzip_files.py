#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2014 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 21 April 2014 13:03 PDT (-0700)
"""


import os
import shutil
import argparse
import ConfigParser

from phyluce.helpers import is_file, is_dir, FullPaths, CreateDir
from phyluce.log import setup_logging

import pdb

def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
        description="""Given a config file and inputs, merge multiple gz files into one"""
    )
    parser.add_argument(
        "--config",
        required=True,
        type=is_file,
        action=FullPaths,
        help="""The path to the config file to use for merging."""
    )
    parser.add_argument(
        "--output",
        required=True,
        action=CreateDir,
        help="""The path to a directory in which to store the output."""
    )
    parser.add_argument(
        "--section",
        type=str,
        help="""The section holding the merge info."""
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
    conf = ConfigParser.ConfigParser()
    conf.optionxform = str
    conf.read(args.config)
    items = conf.items("samples")
    #pdb.set_trace()
    for item in items:
        name, file_names = item
        files = file_names.strip().split(",")
        with open(os.path.join(args.output, name), 'wb') as outfile:
            for infile in sorted(files):
                shutil.copyfileobj(open(infile), outfile)
                log.info("Copied {} to {}".format(
                    os.path.basename(infile),
                    name
                ))
    # end
    text = " Completed {} ".format(my_name)
    log.info(text.center(65, "="))




if __name__ == '__main__':
    main()
