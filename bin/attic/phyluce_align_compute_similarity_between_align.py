#!/usr/bin/env python
# encoding: utf-8
"""
File: compute_similarity_between_align.py
Author: Brant Faircloth

Created by Brant Faircloth on 05 April 2012 15:04 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description:

"""


import os
import re
import sys
import glob
import argparse

from collections import defaultdict
from multiprocessing import cpu_count, Pool

from Bio import AlignIO
from cogent import LoadSeqs
from cogent.phylo import distance
from cogent.evolve.models import GTR

from phyluce.helpers import is_dir, FullPaths

import pdb


def get_args():
    parser = argparse.ArgumentParser(
            description="""Compute summary parameters for alignments"""
        )
    parser.add_argument(
            '--nexus',
            type=is_dir,
            action=FullPaths,
            required=True,
            help='The directory containing the nexus files'
        )
    parser.add_argument(
            '--output',
            type=is_dir,
            action=FullPaths,
            default=None,
            help='The directory to hold trimmed alignments'
        )
    parser.add_argument(
            "--regex",
            action="store_true",
            default=False,
            help="""Trim alignments using regular expressions""",
        )
    parser.add_argument(
            "--multiprocessing",
            action="store_true",
            default=False,
            help="""Help text""",
        )
    parser.add_argument(
            "--keep-alignments",
            dest="keep",
            action="store_true",
            default=False,
            help="""Help text""",
        )
    return parser.parse_args()


def get_files(input_dir):
    return glob.glob(os.path.join(os.path.expanduser(input_dir), '*.nex*'))


def get_aln_starts(aln):
    mn, mx = 0, aln.get_alignment_length()
    for sp in aln:
        for count, base in enumerate(sp.seq):
            if base != '-':
                tf = count
                break
        for count, base in enumerate(sp.seq[::-1]):
            if base != '-':
                tr = count
                break
        tr = len(sp.seq) - tr
        if tf > mn:
            mn = tf
        if mx == 0:
            mx = tr
        elif tr < mx:
            mx = tr
    return aln[:, mn:mx]


def trim_with_regex(aln, forw, rev):
    mn, mx = 0, aln.get_alignment_length()
    for sp in aln:
        f = forw.search(str(sp.seq))
        e = rev.search(str(sp.seq))
        if f and f.end() > mn:
            mn = f.end()
        if e:
            if mx == 0:
                mx = e.start()
            elif e.start() < mx:
                mx = e.start()
    return aln[:, mn:mx]


def get_sequence_divergence(aln, forw=None, rev=None):
    aln = get_aln_starts(aln)
    if forw and rev:
        aln = trim_with_regex(aln, forw, rev)
    diff_gaps = []
    diff_no_gaps = []
    len_included_bases = []
    for column in xrange(aln.get_alignment_length()):
        #pdb.set_trace()
        bases = set([base.upper() for base in aln[:, column]])
        if '?' in bases or 'N' in bases or 'X' in bases:
            pass
        else:
            len_included_bases.append(1)
            if len(bases) > 1:
                diff_gaps.append(1)
            bases.discard('-')
            if len(bases) > 1:
                diff_no_gaps.append(1)
    #pdb.set_trace()
    return aln, sum(len_included_bases), sum(diff_gaps), sum(diff_no_gaps)


def get_pairwise_distance(aln):
    seqs = []
    for sp in aln:
        seqs.extend(['>{}'.format(sp.id), str(sp.seq)])
    paln = LoadSeqs(data=seqs)
    d = distance.EstimateDistances(paln, submodel=GTR())
    d.run(show_progress=False)
    pd = d.getPairwiseDistances()
    return pd.values()


def worker(work):
    f, forw, rev = work
    sys.stdout.write('.')
    sys.stdout.flush()
    diverge, dist = defaultdict(list), defaultdict(list)
    locus = os.path.basename(f)
    sys.stdout.write('.')
    sys.stdout.flush()
    aln = AlignIO.read(f, 'nexus')
    aln, len_included_bases, div, divnogap = get_sequence_divergence(aln, forw, rev)
    diverge[locus] = [
            len_included_bases,
            div,
            divnogap
        ]
    dist[locus] = get_pairwise_distance(aln)
    return [diverge, dist, aln]


def get_sum_diffs(results):
    diff = []
    for r in results:
        diff.append(r[0].values()[0][1])
    return sum(diff)


def get_sum_lengths(results):
    lengths = []
    for r in results:
        lengths.append(r[0].values()[0][0])
    return sum(lengths)


def get_average_distance(results):
    dist = []
    for r in results:
        dist.append(r[1].values()[0][0])
    return float(sum(dist)) / len(dist)


def main():
    args = get_args()
    if args.regex:
        # really aggresive trimming
        #for_regex = re.compile('^[ACGTNacgtn]{0,5}-+[ACGTNacgtn]{,20}-{10,}[ACGTNacgtn]{,20}-{10,}')
        #rev_regex = re.compile('-{10,}[ACGTNacgtn]{,20}-{10,}[ACGTNacgtn]{,20}-+[ACGTNacgtn]{,5}$')
        # less aggressive trimming
        print "Trimming all alignments using regular expressions...\n"
        for_regex = re.compile('^[ACGTNacgtn]{0,10}-+[ACGTNacgtn]{0,5}-{10,}[ACGTNacgtn]{0,5}-{10,}')
        rev_regex = re.compile('-{10,}[ACGTNacgtn]{0,5}-{10,}[ACGTNacgtn]{0,5}-+[ACGTNacgtn]{0,10}$')
    else:
        for_regex, rev_regex = None, None
    # get files and packge for map()
    files = get_files(args.nexus)
    files = [[f, for_regex, rev_regex] for f in files]
    if args.multiprocessing:
        pool = Pool(cpu_count() - 1)
        results = pool.map(worker, files)
    else:
        results = map(worker, files)
    sum_diffs = get_sum_diffs(results)
    total_bases = get_sum_lengths(results)
    average_distance = get_average_distance(results)
    print "\n"
    print "sequence divergence = {0:.2f}%".format(
        float(sum_diffs) / total_bases * 100
        )
    print "average pairwise distance = {0:.4f}".format(average_distance)
    divergewriter = open(os.path.join(args.output, 'sequence_divergence.csv'), 'w')
    divergewriter.write('locus,len,diff(gaps),diff(nogaps)\n')
    distwriter = open(os.path.join(args.output, 'pairwise_distance.csv'), 'w')
    distwriter.write('locus,pairwise_distance\n')
    for r in results:
        diverge, dist, aln = r
        for locus, values in diverge.iteritems():
            divergewriter.write("{},{},{},{}\n".format(locus, values[0], values[1], values[2]))
        for locus, values in dist.iteritems():
            distwriter.write("{},{}\n".format(locus, values[0]))
            if args.keep:
                aln_out = open(os.path.join(args.output, locus), 'w')
                AlignIO.write(aln, aln_out, 'nexus')
    for i in [divergewriter, distwriter]:
        i.close()


if __name__ == '__main__':
    main()
