#!/usr/bin/env python
# encoding: utf-8
"""
File: unnamed_file.py
Author: Brant Faircloth

Created by Brant Faircloth on 31 August 2013 13:08 PDT (-0700)
Copyright (c) 2013 Brant C. Faircloth. All rights reserved.

Description:

"""

import os
import sys
import glob
import argparse
import subprocess
import multiprocessing
from phyluce.helpers import FullPaths, is_dir, is_file
from Bio import SeqIO

import pdb

def get_args():
    """Run SATé alignments in parallel"""
    parser = argparse.ArgumentParser(
        description="""parallel_sate.py"""
    )
    parser.add_argument(
        "--input",
        required=True,
        action=FullPaths,
        type=is_dir,
        help="""The input directory containing fasta files"""
    )
    parser.add_argument(
        "--output",
        required=True,
        action=FullPaths,
        type=is_dir,
        help="""The output directory to hold alignments"""
    )
    parser.add_argument(
        "--sate-conf",
        required=True,
        type=is_file,
        help="""The path to the SATé config file"""
    )
    parser.add_argument(
        "--cores",
        type=int,
        default=1,
        help="""The number of compute cores to use"""
    )
    return parser.parse_args()

def get_files(input_dir):
    extensions = ('.fasta', '.fas', '.fsa')
    files = []
    for ext in extensions:
        files.extend(glob.glob(os.path.join(os.path.expanduser(input_dir), '*{}*'.format(ext))))
    # ensure we collapse duplicate filenames
    unique_files = list(set(files))
    return_files = []
    # drop fasta files with < 2 taxa
    for f in unique_files:
        with open(f, 'rU') as infile:
            record = list(SeqIO.parse(infile, 'fasta'))
            if len(record) >= 4:
                return_files.append(f)
            else:
                print "Dropping {} because < 4 taxa".format(os.path.basename(f))
    return return_files

def run_sate(work):
    #pdb.set_trace()
    file, output, config = work
    # get locus name from path
    jobname = os.path.basename(file).split('.')[0]
    cmd = [
    'run_sate.py',
    '--input',
    file,
    '--output-directory',
    output,
    '--job',
    jobname,
    '--missing',
    'Ambiguous',
    config
    ]
    proc = subprocess.Popen(cmd,
        stderr=subprocess.PIPE,
        stdout=subprocess.PIPE
    )
    stdout, stderr = proc.communicate()
    #keep = set()
    #files = glob.glob("{}*".format(jobname))
    #for i in files:
    #    for ext in [".err.txt",".out.txt", ".aln"]:
    #        if i.endswith(ext):
    #            keep.add(i)
    #remove = set(files).difference(keep)
    #for r in remove:
    #    os.remove(os.path.join(output, r))
    if not stderr:
        sys.stdout.write(".")
        sys.stdout.flush()
        return None
    else:
        return jobname

def main():
    args = get_args()
    fastas = get_files(args.input)
    work = [[f, args.output, args.sate_conf] for f in fastas]
    if args.cores > 1:
        assert args.cores <= multiprocessing.cpu_count(), "You've specified more cores than you have"
        pool = multiprocessing.Pool(args.cores)
        results = pool.map(run_sate, work)
    else:
        results = map(run_sate, work)
    if set(results) != set([None]):
        print "There may be a problem in:"
        for i in results:
            if i != None:
                print i

if __name__ == '__main__':
    main()
