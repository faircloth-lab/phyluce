#!/usr/bin/env python
# encoding: utf-8
"""
File: randomly_sample_and_concatenate.py
Author: Brant Faircloth

Created by Brant Faircloth on 26 April 2012 15:04 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: Given an input directory of nexus files, randomly sample
some, concatenate them, and return a conatenated nexus file as output.

"""

import os
import sys
import glob
import numpy
import argparse
from Bio.Nexus import Nexus
from phyluce.helpers import is_dir, FullPaths

import pdb


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Randomly sample a directory of nexus files, concatenate, and output""")
    parser.add_argument(
            "nexus",
            action=FullPaths,
            type=is_dir,
            help="""The input directory of nexus files"""
        )
    parser.add_argument(
            "output",
            action=FullPaths,
            type=is_dir,
            help="""The output directory of contenated files"""
        )
    parser.add_argument(
            "--sample-size",
            dest='sample_size',
            type=int,
            default=50,
            help="""The number of loci to sample"""
        )
    parser.add_argument(
            "--samples",
            type=int,
            default=1,
            help="""The number of samples to take""",
        )
    return parser.parse_args()


def main():
    args = get_args()
    #pdb.set_trace()
    # get filenames in directory and convert to array
    files = numpy.array(glob.glob(os.path.join(args.nexus, '*.nex*')))
    # make sure we have enough
    assert len(files) >= args.sample_size, "Sample size must be < number(files)"
    print "Running"
    for i in xrange(args.samples):
        sys.stdout.write('.')
        sys.stdout.flush()
        # get list of random numbers
        sample = numpy.random.random_integers(0, len(files) - 1, args.sample_size)
        # reindex filenames by random selections
        random_files = sorted(files[sample].tolist())
        # concatenate and output
        files_to_combine = [(f, Nexus.Nexus(f)) for f in random_files]
        combined = Nexus.combine(files_to_combine)
        align_name = "random-sample-{}-{}-loci.nex".format(i, args.sample_size)
        # open metadata file
        meta_name = 'META-random-sample-{}-{}-loci.txt'.format(i, args.sample_size)
        meta = open(
                os.path.join(args.output, meta_name), 'w'
            )
        meta.write('{}'.format('\n'.join(random_files)))
        meta.close()
        combined.write_nexus_data(filename=open(
                os.path.join(args.output, align_name), 'w')
            )
    sys.stdout.write("Done")

if __name__ == '__main__':
    main()
