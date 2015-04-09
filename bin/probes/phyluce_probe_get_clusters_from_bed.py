#!/usr/bin/env python
# encoding: utf-8
"""
File: get_clusters_from_bed.py
Author: Brant Faircloth

Created by Brant Faircloth on 17 July 2012 16:07 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: Given an input BED file, find clusters of overlapping
entries.

"""

import os
import glob
import json
import argparse
from collections import Counter, defaultdict
from bx.intervals.cluster import ClusterTree
from phyluce.helpers import FullPaths

import pdb


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Program description""")
    parser.add_argument(
            "bed",
            action=FullPaths,
            help="""The BED directory you want to search for clusters"""
        )
    parser.add_argument(
            "outdir",
            action=FullPaths,
            help="""The output directory to store results"""
        )
    parser.add_argument(
            "--files",
            nargs='+',
            default=None,
            help="""Specific files in the directory you want to process""",
        )
    return parser.parse_args()


def get_dict_from_name(name):
    return dict([n.split(':') for n in name.split('|')])


def check_for_dupe_hits(counter):
    return {name:count for name, count in counter.iteritems() if count > 1}


def main():
    args = get_args()
    if args.files:
        files = [f for f in glob.glob(os.path.join(args.bed, '*.bed')) if f in args.files]
        assert len(files) == len(args.files), "You have specified files that are not in {0}".format(args.bed)
    else:
        files = glob.glob(os.path.join(args.bed, '*.bed'))
    #pdb.set_trace()
    for f in files:
        # setup output files
        fname = os.path.splitext(os.path.basename(f))[0]
        print "Processing: {0}".format(fname)
        outfname = "{0}.json".format(fname)
        outfile = open(os.path.join(args.outdir, outfname), 'w')
        # outdata
        outdata = {}
        # setup our cluster tree, a name dict, and a counter
        cluster_trees = defaultdict(lambda: ClusterTree(500, 2))
        name_map = {}
        name_counter = Counter()
        for line in open(f, 'rU'):
            if not line.startswith('track'):
                chromo, start, end, name = line.split("\t")[:4]
                start, end = int(start), int(end)
                # convert the BED name to a dict, indexed by unique db id pkey
                temp_dict = get_dict_from_name(name)
                key = int(temp_dict['probes-id'])
                # create a counter of names, so we can check for dupe hits
                name_counter.update([key])
                name_map[key] = temp_dict
                cluster_trees[chromo].insert(start, end, key)
        duplicate_hits = check_for_dupe_hits(name_counter)
        if duplicate_hits:
            outdata['Duplicates'] = True
        else:
            outdata['Duplicates'] = False
        outdata['Overlaps'] = defaultdict(list)
        for chromo in cluster_trees:
            overlaps = cluster_trees[chromo].getregions()
            if overlaps:
                for span in overlaps:
                    # get distinct list of *loci* (not probes) hit for a given
                    # region.
                    loci = list(set([int(name_map[probe]['probes-locus']) for probe in span[2]]))
                    if len(span[2]) > 1:
                        outdata['Overlaps'][chromo].append({
                                'start': span[0],
                                'end': span[1],
                                'loci': loci,
                                'probes': span[2]
                            })
        json.dump(outdata, outfile, indent=2)


if __name__ == '__main__':
    main()


