#!/usr/bin/env python
# encoding: utf-8
"""
File: run_raxml_genetrees.py
Author: Brant Faircloth

Created by Brant Faircloth on 13 September 2012 18:09 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description:

"""

import os
import re
import sys
import glob
import random
import argparse
import subprocess
import multiprocessing
from phyluce.helpers import is_dir, FullPaths
from phyluce.third_party import which

import pdb

def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Program description""")
    parser.add_argument(
            "--input",
            required=True,
            type=is_dir,
            action=FullPaths,
            help="""The input directory containing alignments in phylip format"""
        )
    parser.add_argument(
            "--output",
            required=True,
            type=is_dir,
            action=FullPaths,
            help="""The output directory to hold alignments"""
        )
    parser.add_argument(
            "--outgroup",
            type=str,
            help="""The outgroup to use"""
        )
    parser.add_argument(
            "--threads",
            type=int,
            default=1,
            help="""The number of RAxML threads to run (best to determine empirically)"""
        )
    parser.add_argument(
            "--cores",
            type=int,
            default=1,
            help="""The number of concurrent RAxML jobs to run"""
        )
    parser.add_argument(
            "--quiet",
            action="store_true",
            default=False,
            help="""Suppress the CPU usage question""",
        )
    return parser.parse_args()


def get_many_thread_raxml_cmd(threads, outgroup, alignment, outputdir):
    cmd = [
            which("raxmlHPC-PTHREADS-SSE3")[0],
            "-T",
            str(threads),
            "-o",
            outgroup,
            "-m",
            "GTRCAT",
            "-n",
            "best",
            "-s",
            alignment,
            "-N",
            '10',
            "-p",
            str(random.randint(0, 1000000)),
            "-w",
            outputdir
        ]
    return cmd


def get_one_thread_raxml_cmd(threads, outgroup, alignment, outputdir):
    cmd = [
            which("raxmlHPC-SSE3")[0],
            "-o",
            outgroup,
            "-m",
            "GTRCAT",
            "-n",
            "best",
            "-s",
            alignment,
            "-N",
            '10',
            "-p",
            str(random.randint(0, 1000000)),
            "-w",
            outputdir
        ]
    return cmd


def run_raxml(work):
    threads, output, outgroup, time, patterns, alignment = work
    # get the alignment name
    dirname = os.path.splitext(os.path.basename(alignment))[0]
    # make a directory for the alignment; raxml needs trailing slash
    outputdir = os.path.join(output, dirname) + "/"
    os.makedirs(outputdir)
    # run raxml saving the output to a directory in args.output named
    # after the alignment:
    #
    # raxmlHPC-PTHREADS-SSE3 -T {threads} -o {outgroup} -m GTRCAT -n best -s {alignment} -N 10 -p {random #} -w {outputdir}
    if threads == 1:
        cmd = get_one_thread_raxml_cmd(threads, outgroup, alignment, outputdir)
    else:
        cmd = get_many_thread_raxml_cmd(threads, outgroup, alignment, outputdir)
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    #pdb.set_trace()
    # get run time for number of patterns
    seconds = time.search(stdout).groups()[0]
    sites = patterns.search(stdout).groups()[0]
    sys.stdout.write("name={0},seconds={1},sites={2}\n".format(dirname, seconds, sites))
    sys.stdout.flush()
    return outputdir


def main():
    args = get_args()
    # get the number of jobs as available_procs / threads for raxml
    jobs = args.cores / args.threads
    question = "The total number of cores in use is {0}. This will run \n" + \
        "{1} concurrent jobs of {2} thread(s). Is this correct [Y/n]? "
    if args.quiet:
        correct_jobs = "Y"
    else:
        correct_jobs = raw_input(question.format(
            args.cores,
            jobs,
            args.threads)
        )
    if correct_jobs == "Y":
        assert jobs < multiprocessing.cpu_count(), "The total number of jobs * threads is greather than the available CPUs"
        pool = multiprocessing.Pool(jobs)
        time = re.compile("Overall\sexecution\stime:\s(\d+\.\d+)\ssecs")
        patterns = re.compile("Alignment\sPatterns:\s(\d+)")
        alignments = glob.glob(os.path.join(args.input, '*.phylip'))
        work = [[args.threads, args.output, args.outgroup, time, patterns, alignment] for alignment in alignments]
        trees = pool.map(run_raxml, work)
        output = open(os.path.join(args.output, "all-best-trees.tre"), 'w')
        for treedir in trees:
            best_tree = open(os.path.join(treedir, "RAxML_bestTree.best"), 'rb').read()
            output.write(best_tree)
        output.close()

if __name__ == '__main__':
    main()
