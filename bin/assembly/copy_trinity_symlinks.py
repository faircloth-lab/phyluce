#!/usr/bin/env python
# encoding: utf-8
"""
File: copy_trinity_symlinks.py
Author: Brant Faircloth

Created by Brant Faircloth on 08 June 2013 17:06 PDT (-0700)
Copyright (c) 2013 Brant C. Faircloth. All rights reserved.

Description:

"""

import os
import sys
import argparse
import ConfigParser
from phyluce.helpers import FullPaths

import pdb


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
        description="""Copy symlinks from trinity assemblies to a new folder"""
    )
    parser.add_argument(
        "--assembly-symlinks",
        dest="symlinks",
        required=True,
        action=FullPaths,
        help="""The location of the trinity symlinks (trinity-assemblies/contigs)"""
    )
    parser.add_argument(
        "--conf",
        help="""The configuration file to use"""
    )
    parser.add_argument(
        "--output",
        action=FullPaths,
        help="""The output folder"""
    )
    return parser.parse_args()


def main():
    args = get_args()
    config = ConfigParser.RawConfigParser(allow_no_value=True)
    config.read(args.conf)
    sections = config.sections()
    links = ()
    # check for files and create links first so we don't make a
    # half-empty folder of symlinks
    for section in sections:
        if section in ['outgroups', 'taxa']:
            taxa = config.items(section)
            for taxon in taxa:
                link_name = "{}.contigs.fasta".format(taxon[0])
                link_path = os.path.join(args.symlinks, link_name)
                if os.path.islink(link_path):
                    link = os.readlink(link_path)
                    dst = os.path.join(args.output, link_name)
                    links += ((link, dst),)
                else:
                    raise OSError("The path {} does not exist".format(link_path))
    # Now that we've found the links and validated they exist,
    # make the output directory
    try:
        os.makedirs(args.output)
    except OSError:
        answer = raw_input("Output directory exists.  Overwrite [Y/n]? ")
        if answer != 'Y':
            sys.exit()
    # make the links
    for src, dst in links:
        try:
            os.symlink(src, dst)
        except OSError:
            raise OSError("The path {} already exists".format(dst))


if __name__ == '__main__':
    main()
