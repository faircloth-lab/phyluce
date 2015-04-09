#!/usr/bin/env python
# encoding: utf-8
"""
File: thin_mrbayes_runs.py
Author: Brant Faircloth

Created by Brant Faircloth on 29 March 2012 09:03 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: 

"""

import os
import sys
import glob
import shutil
import argparse
from phyluce.helpers import FullPaths, is_dir

import pdb

def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Thin a folder of mrbayes output""")
    parser.add_argument(
            "input",
            action=FullPaths,
            type=is_dir,
            help="""Input directory"""
        )
    parser.add_argument(
            "output",
            action=FullPaths,
            help="""Output directory"""
        )
    parser.add_argument(
            "--thin",
            type=int,
            default=100,
            help="""Thinning factor""",
        )
    return parser.parse_args()

def thin_p_files(input, output, thin = 100):
    for line in open(input, 'rU'):
        if not line.split('\t')[0].isdigit():
            output.write(line)
        else:
            if int(line.split('\t')[0]) == 1:
                output.write(line)
            elif int(line.split('\t')[0]) % thin == 0:
                output.write(line)

def thin_t_files(input, output, thin = 100):
    for line in open(input, 'rU'):
        if not line.startswith('   tree rep.'):
            output.write(line)
        else:
            #pdb.set_trace()
            if int(line.split('=')[0].strip().split('.')[1]) == 1:
                output.write(line)
            elif int(line.split('=')[0].strip().split('.')[1]) % thin == 0:
                output.write(line)

def main():
    args = get_args()
    files = [f for f in glob.glob(os.path.join(args.input, '*')) if os.path.splitext(f)[1] in ['.t', '.p']]
    # check if outputdir
    try:
        assert os.path.exists(args.output)
    except:
        inp = raw_input('Output directory does not exist.  Create [Y/n]: ')
        if inp == 'Y':
            os.makedirs(args.output)
        else:
            print "Exiting"
            sys.exit()
    nexus = glob.glob(os.path.join(args.input, '*.nex'))
    assert len(nexus) == 1, "There is more than one nexus file"
    output_name = os.path.basename(nexus[0])
    output_file = os.path.join(args.output, output_name)
    shutil.copyfile(nexus[0], output_file)
    for input in files:
        output_name = os.path.basename(input)
        output_file = os.path.join(args.output, output_name)
        output = open(output_file, 'w')
        if os.path.splitext(input)[1] == '.p':
            thin_p_files(input, output, args.thin)
        elif os.path.splitext(input)[1] == '.t':
            thin_t_files(input, output, args. thin)
        output.close()


if __name__ == '__main__':
    main()