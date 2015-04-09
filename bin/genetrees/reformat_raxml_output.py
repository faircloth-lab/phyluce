#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2014 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 12 April 2014 15:22 PDT (-0700)
"""

import os
import glob
import argparse
from collections import defaultdict
from phyluce.helpers import is_dir, FullPaths

#import pdb

def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
        description="""Program description"""
    )
    parser.add_argument(
        "--genetrees",
        action=FullPaths,
        type=is_dir,
        help="""The path to the directory containing RAxML best trees"""
    )
    parser.add_argument(
        "--bootreps",
        action=FullPaths,
        type=is_dir,
        help="""The path to the directory containing RAxML bootreps"""
    )
    parser.add_argument(
        "--output",
        required=True,
        action=FullPaths,
        help="""The output file name"""
    )
    return parser.parse_args()


def main():
    args = get_args()
    if args.genetrees:
        with open("{0}.genetrees".format(args.output), 'w') as outf:
            for file in glob.glob(os.path.join(args.genetrees, "RAxML_bestTree.*")):
                with open(file, "rU") as infile:
                    tree = infile.read().strip()
                outf.write("{0}\n".format(tree))
    elif args.bootreps:
        tree_dict = defaultdict(list)
        for file in glob.glob(os.path.join(args.bootreps, "RAxML_bootstrap.*")):
            #pdb.set_trace()
            with open(file, "rU") as infile:
                for cnt, line in enumerate(infile):
                    tree_dict[cnt].append(line.strip())
        with open("{0}.bootreps".format(args.output), 'w') as outf:
            for cnt, trees in tree_dict.iteritems():
                for tree in trees:
                    outf.write("{0}\t{1}\n".format(cnt, tree))



if __name__ == '__main__':
    main()
