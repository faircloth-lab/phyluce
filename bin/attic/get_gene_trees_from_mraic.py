#!/usr/bin/env python
# encoding: utf-8
"""
File: get_gene_trees_from_mraic.py
Author: Brant Faircloth

Created by Brant Faircloth on 18 April 2012 22:04 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: Used to generate a cloudforest-like *.tre file
from run_mraic.py (often used to complete "failed" jobs)

"""

import os
import sys
import glob
import argparse

from phyluce.helpers import is_dir, FullPaths

import pdb

def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Generate a cloudforest-like genetree file from run_mraic output""")
    parser.add_argument(
            "mraic",
            action=FullPaths,
            type=is_dir,
            help="""The run_mraic output directory"""
        )
    parser.add_argument(
            "output",
            type=argparse.FileType('w'),
            help="""The output file in which to store genetrees"""
        )
    parser.add_argument(
            "--model-type",
            dest="model",
            type=str,
            choices=["AICc", "AIC", "BIC"],
            default="AICc",
            help="""The model type to parse from the results folder""",
        )
    return parser.parse_args()

def make_tree_name(name, model):
    return "chrm={},model={}".format(name, model)

def main():
    args = get_args()
    for f in glob.glob(os.path.join(args.mraic, "*.{}-*".format(args.model))):
        fs = os.path.basename(f).split('.')
        name = fs[0]
        model = fs[2].strip("{}-".format(args.model))
        tree = open(f, 'rU').read().strip()
        tree = "tree '" + make_tree_name(name, model) + "' = [&U] " + tree
        args.output.write('null\t"{}"\n'.format(tree))
    args.output.close()

if __name__ == '__main__':
    main()
