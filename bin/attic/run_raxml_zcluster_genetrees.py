#!/usr/bin/env python
# encoding: utf-8
"""
File: run_raxml_genetrees_on_zcluster.py
Author: Brant Faircloth

Created by Brant Faircloth on 07 March 2013 12:03 PST (-0800)
Copyright (c) 2013 Brant C. Faircloth. All rights reserved.

Description: 

"""

import os
import glob
#import time
import random
import argparse
import subprocess
from phyluce.helpers import FullPaths, is_dir, get_file_extensions


import pdb


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Generate raxml genetrees on zcluster (qsub) system""")
    parser.add_argument(
            "input",
            action=FullPaths,
            type=is_dir,
            help="""The path holding the alignments"""
        )
    parser.add_argument(
            "--output",
            default="raxml-output",
            help="""The directory to hold output files"""
        )
    parser.add_argument(
            "--control",
            default="submit-scripts",
            help="""The directory to hold control files"""
        )
    parser.add_argument(
            "--outgroup",
            type=str,
            help="""Outgroup name""",
        )
    parser.add_argument(
            "--input-format",
            dest="input_format",
            choices=['fasta', 'nexus', 'phylip', 'clustal', 'emboss', 'stockholm'],
            default='phylip',
            help="""The input alignment format""",
        )
    return parser.parse_args()


def get_files(input_dir, input_format):
    alignments = []
    for ftype in get_file_extensions(input_format):
        alignments.extend(glob.glob(os.path.join(input_dir, "*{}".format(ftype))))
    return alignments


def main():
    args = get_args()
    cwd = os.getcwd()
    newdir = os.path.join(cwd, "raxml-output")
    controldir = os.path.join(cwd, "submit-scripts")
    # make a directory to hold the output
    os.mkdir(newdir)
    os.mkdir(controldir)
    # change to the controldir
    os.chdir(controldir)
    # iterate over files, get file name
    for f in get_files(args.input, args.input_format):
        name = os.path.splitext(os.path.split(f)[1])[0]
        script = name + "_control.sh"
        random_integer = random.randint(0, 1000000)
        cli = "#!/bin/bash\ncd {0}\n/usr/local/raxml/latest/raxmlHPC-SSE3 -o {1} -m GTRCAT -n {2}.best -s {3} -N 10 -p {4} -w {5}".format(
            controldir,
            args.outgroup,
            name,
            f,
            random_integer,
            newdir
            )
        control = open(script, 'w')
        control.write(cli)
        control.close()
        # run SATe
        qsub = [
                'qsub',
                '-q',
                "rcc-30d",
                script
            ]
        stderr, stdout = subprocess.Popen(
                qsub,
                stderr=subprocess.PIPE,
                stdout=subprocess.PIPE
            ).communicate()
        print stdout, stderr
        break


if __name__ == '__main__':
    main()
