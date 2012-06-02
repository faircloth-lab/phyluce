#!/usr/bin/env python
# encoding: utf-8
"""
File: concatenate_gene_trees_from_cloudforest.py
Author: Brant Faircloth

Created by Brant Faircloth on 17 April 2012 22:04 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: 

"""

import os
import sys
import glob
import shutil
import argparse

from phyluce.helpers import is_dir, FullPaths

import pdb


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Concatenate genetrees output from cloudforest""")
    parser.add_argument(
            "genetrees",
            action=FullPaths,
            type=is_dir,
            help="""The directory containing gene trees to concatenate"""
        )
    parser.add_argument(
            "output",
            type=argparse.FileType('w'),
            default=sys.stdout,
            help="""The directory containing gene trees to concatenate"""
        )
    parser.add_argument(
            "--alignments",
            action=FullPaths,
            type=is_dir,
            default=False,
            help="""The directory containing alignments to ensure each has a gene tree""",
        )
    parser.add_argument(
            "--models",
            type=argparse.FileType('w'),
            default=False,
            help="""The output file to save substitution models for each locus""",
        )
    return parser.parse_args()


def get_alignments(alignments):
    if alignments:
        alns = set([os.path.splitext(os.path.basename(f))[0] for f in glob.glob(os.path.join(alignments, "*"))])
    else:
        alns = False
    return alns


def get_and_write_genetrees(genetrees, output):
    loci = []
    models = []
    for f in glob.glob(os.path.join(genetrees, "*")):
        for line in open(f, 'rU'):
            output.write(line)
            # parse genetree lines for locus name
            name, model = line.split(' ')[1].strip("'").split(',')
            locus = name.split('=')[1]
            subs = model.split('=')[1]
            loci.append(locus)
            models.append([locus, subs])
    assert len(loci) == len(set(loci)), "There are duplicate loci"
    assert len(models) == len(set([m[0] for m in models])), "There duplicate models"
    assert len(models) == len(loci), "There different numbers of models and loci"
    return loci, models


def write_models(models, output):
    for model in models:
        output.write("{}\tAICc-{}\n".format(model[0], model[1]))

def copy_over_missing_loci(alignments, diff):
    os.makedirs('missing-loci')
    dst = os.path.abspath('missing-loci')
    for f in glob.glob(os.path.join(alignments, "*")):
        name = os.path.splitext(os.path.basename(f))[0]
        if name in diff:
            shutil.copyfile(f, os.path.join(dst, os.path.basename(f)))

def main():
    args = get_args()
    alns = get_alignments(args.alignments)
    loci, models = get_and_write_genetrees(args.genetrees, args.output)
    if args.models:
        write_models(models, args.models)
    diff = alns.difference(set(loci))
    if diff:
        copy_over_missing_loci(args.alignments, diff)


if __name__ == '__main__':
    main()
