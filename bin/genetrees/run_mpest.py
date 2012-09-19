#!/usr/bin/env python
# encoding: utf-8

"""
run_mpest.py

Created by Brant Faircloth on 16 September 2010 09:37 PDT (-0700).
Copyright (c) 2010 Brant C. Faircloth. All rights reserved.

PURPOSE:  Run an input file through our slightly customized version of mpest
either using a single core or multiple cores, and feeding mpest a random
(integer) seed drawn from a uniform distribution.

USAGE:  python run_mpest.py control-file-name.control output-file-name.output \
            --iterations 1000 --cores 7

Modified MP-EST is at:  https://github.com/faircloth-lab/mp-est

"""


import os
import re
import sys
import time
import shutil
import random
import argparse
import tempfile
#import dendropy
import subprocess
from operator import itemgetter
from multiprocessing import Pool
from collections import defaultdict
from phyluce.helpers import is_file, FullPaths
from ete2 import Tree

import pdb


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Run MP-EST using multiple compute cores""")
    parser.add_argument(
            "control",
            type=is_file,
            action=FullPaths,
            help="""The MP-EST control file"""
        )
    parser.add_argument(
            "output",
            action=FullPaths,
            help="""The MP-EST output file (to hold iterations"""
        )
    parser.add_argument(
            "root",
            type=str,
            help="""The nodename on which to root the tree"""
        )
    parser.add_argument(
            "--iterations",
            type=int,
            default=1,
            help="""The number of iterations to run""",
        )
    parser.add_argument(
            "--cores",
            type=int,
            default=1,
            help="""The number of compute cores to use""",
        )
    parser.add_argument(
            "--bootreps",
            action="store_true",
            default=False,
            help="""If processing bootreps""",
        )
    parser.add_argument(
            "--bootrep-num",
            dest='bootrep_num',
            type=int,
            default=10,
            help="""The number of bootreps to run""",
        )
    parser.add_argument(
            "--raxml",
            action="store_true",
            default=False,
            help="""If processing raxml input (versus cloudforest)""",
        )
    return parser.parse_args()


def create_control_file(taxa, genes, temp_bootrep_file):
    temp_fd, temp_out = tempfile.mkstemp(suffix='.mpest-control')
    taxa_string = ['{0}\t1\t{0}'.format(name) for name in taxa]
    data = {
            'filename': temp_bootrep_file,
            'tree_count': genes,
            'numb_taxa': len(taxa),
            'taxa_details': '\n'.join(taxa_string)
        }
    template = "%(filename)s\n0\n%(tree_count)s %(numb_taxa)s\n%(taxa_details)s\n0\n" % data
    os.write(temp_fd, template)
    os.close(temp_fd)
    return temp_out


def run_mpest(work):
    taxa, genes, temp_bootrep_file = work
    #pdb.set_trace()
    # generate control file
    control = create_control_file(taxa, genes, temp_bootrep_file)
    # create output file
    temp_fd, temp_out = tempfile.mkstemp(suffix='.mpest-out')
    os.close(temp_fd)
    # get seed
    seed = random.randint(1, 1000000)
    # setup mpest command line
    cmd = ['mpest', control, str(seed), temp_out]
    # capture output to hide it from CLI
    out, err = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        ).communicate(None)
    regex = re.compile("tree\smpest\s\[(.*)\]")
    for line in reversed(open(temp_out, 'rU').readlines()):
        # we basically need to pull the last line of the file
        if 'tree mpest' in line:
            pre, post = line.strip('\n').split('=')
            post = post.strip()
            ll_result = regex.search(pre)
            ll = float(ll_result.groups()[0])
            result = post.rstrip(';')
            break
    #for file in [control, temp_out, temp_bootrep_file]:
    #    os.remove(file)
    sys.stdout.write('.')
    sys.stdout.flush()
    return result


def parse_bootrep_file(fname, root, bootrep_num):
    bootrep_temp_files = []
    #bootreps = defaultdict(dendropy.TreeList)
    bootreps = defaultdict(list)
    sys.stdout.write("Parsing bootrep file")
    sys.stdout.flush()
    printrep = 1
    for line in open(fname, 'rU'):
        repnum, tree_string = line.strip().split('\t')
        # clean up input
        repnum = int(repnum)
        if repnum > printrep:
            sys.stdout.write('.')
            sys.stdout.flush()
            printrep = repnum
        tree_string = tree_string.strip('"')
        if repnum <= bootrep_num:
            tree = Tree(tree_string)
            tree.set_outgroup(root)
            bootreps[repnum].append(tree)
        else:
            break
    genes = list(set([len(trees) for trees in bootreps.values()]))
    assert len(genes) == 1
    for repnum, trees in bootreps.iteritems():
        temp_fd, temp_out = tempfile.mkstemp(prefix='{}-'.format(repnum), suffix='.mpest-bootrep')
        for tree in trees:
            os.write(temp_fd, tree.write(format=5) + "\n")
        os.close(temp_fd)
        bootrep_temp_files.append(temp_out)
    return genes[0], bootrep_temp_files


def parse_genetree_file(fname, root, raxml=False):
    genetree_temp_files = []
    genetrees = defaultdict(list)
    for line in open(fname, 'rU'):
        if raxml:
            tree_string = line.strip()
        else:
            repnum, tree_string = line.strip().split('\t')
        tree = Tree(tree_string)
        tree.set_outgroup(root)
        genetrees[0].append(tree)
    genes = len(genetrees[0])
    for repnum, trees in genetrees.iteritems():
        temp_fd, temp_out = tempfile.mkstemp(prefix='{}-'.format(repnum), suffix='.mpest-bootrep')
        for tree in trees:
            os.write(temp_fd, tree.write(format=5) + "\n")
        os.close(temp_fd)
        genetree_temp_files.append(temp_out)
    return genes, genetree_temp_files


def get_taxa_for_one_alignment(fname, raxml=False):
    line = open(fname, 'rU').readline()
    if raxml:
        tree_string = line.strip()
    else:
        repnum, tree_string = line.strip().split('\t')
    tree_string = tree_string.strip('"')
    tree = Tree(tree_string)
    taxa = tuple(tree.get_leaf_names())
    return taxa


def main():
    start_time = time.time()
    print 'Started: ', time.strftime("%a %b %d, %Y  %H:%M:%S", time.localtime(start_time))
    args = get_args()
    if not args.bootreps:
        # get taxa for a given alignment
        taxa = get_taxa_for_one_alignment(args.control, args.raxml)
        # parse genetrees from output file
        genes, temp_genetree_file = parse_genetree_file(args.control, args.root, args.raxml)
        # create local control file
        local_control = create_control_file(taxa, genes, temp_genetree_file[0])
        shutil.move(local_control, os.path.join(os.getcwd(), 'genetree-local-control-file.control'))
        work = [[taxa, genes, temp_genetree_file[0]]]
        sys.stdout.write("\nRunning MP-EST against {0} genetrees".format(genes))
        sys.stdout.flush()
        tree = map(run_mpest, work)
        outp = open(args.output, 'w')
        pdb.set_trace()
        outp.write(tree[0])
        outp.write(';\n')
        outp.close()
    elif args.bootreps:
        # get taxa for a given alignment
        taxa = get_taxa_for_one_alignment(args.control)
        # parse bootreps from output file
        genes, temp_bootrep_files = parse_bootrep_file(args.control, args.root, args.bootrep_num)
        # create local control file
        local_control = create_control_file(taxa, genes, temp_bootrep_files[0])
        shutil.move(local_control, os.path.join(os.getcwd(), 'local-control-file.control'))
        # cram taxa, gene count, and bootrep files into a list of size num(iterations)
        work = [[taxa, genes, bootrep] for bootrep in temp_bootrep_files]
        sys.stdout.write("\nRunning {0} bootreps of {1} iterations".format(args.bootrep_num, args.iterations))
        sys.stdout.flush()
        if args.cores == 1:
            trees = map(run_mpest, work)
        else:
            pool = Pool(args.cores)
            trees = pool.map(run_mpest, work)
        outp = open(args.output, 'w')
        outp.write(';\n'.join(trees))
        outp.write(';')
        outp.close()
    end_time = time.time()
    print 'Ended: ', time.strftime("%a %b %d, %Y  %H:%M:%S", time.localtime(end_time))
    print 'Time for execution: ', (end_time - start_time) / 60, 'minutes'


if __name__ == '__main__':
    main()
