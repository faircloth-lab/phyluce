#!/usr/bin/env python
# encoding: utf-8
"""
File: get_clusters_from_taxa.py
Author: Brant Faircloth

Created by Brant Faircloth on 19 July 2012 18:07 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: 

"""

import os
import json
import glob
import sqlite3
import argparse
import itertools
from collections import defaultdict, Counter
from phyluce.helpers import FullPaths

import pdb


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Program description""")
    parser.add_argument(
            "clusters",
            action=FullPaths,
            help="""The directory containing cluster files"""
        )
    parser.add_argument(
            "db",
            action=FullPaths,
            help="""The database to update"""
        )
    parser.add_argument(
            "--taxa",
            nargs='+',
            default=None,
            help="""The taxon overlaps to use""",
        )
    return parser.parse_args()


def get_compare_files(args):
    """Return the list of files"""
    compare_files = {}
    for i in glob.glob(os.path.join(args.clusters, '*.json')):
        name = os.path.splitext(os.path.basename(i))[0]
        if args.taxa and name in args.taxa:
            compare_files[name] = i
        else:
            compare_files[name] = i
    return compare_files


def get_overlapping_data_by_taxon(compare_files):
    """Read overlap data from json into a dict"""
    overlaps = {}
    for k, v in compare_files.iteritems():
        print "Working on {0}...".format(k)
        compare_set = set([])
        cdata = json.load(open(v))
        for chromo, overlap in cdata['Overlaps'].iteritems():
            for hit in overlap:
                compare_set.add(frozenset(hit['probes']))
        overlaps[k] = compare_set
    return overlaps


def get_intersection_of_taxa(overlaps, all_overlaps, cluster=3):
    """Given an input list of overlaps and perfect overlaps, add probe clusters
    to perfect overlaps where the cluster exists in `cluster` number of taxa. Then
    go through final list and merge subsets/supersets of loci."""
    for combo in itertools.combinations(overlaps.keys(), cluster):
        sl = [overlaps[i] for i in combo]
        trio = set.intersection(*sl)
        print combo, len(trio)
        all_overlaps = all_overlaps.union(trio)
    nesteds, skips = set([]), set([])
    # determine supersets/subsets of loci
    for combo in itertools.combinations(all_overlaps, 2):
        if combo[0].issubset(combo[1]) or combo[0].issuperset(combo[1]):
            nesteds.add(combo[0].union(combo[1]))
            skips.add(combo[0])
            skips.add(combo[1])
    # remove what we're merging
    for skip in skips:
        all_overlaps.remove(skip)
    # add back what we merged
    for nested in nesteds:
        all_overlaps.add(nested)
    return all_overlaps


def check_probe_counts_in_clusters(all_overlaps):
    """Ensure that only one probe exists in each cluster"""
    probes = []
    for s in all_overlaps:
        probes.extend(list(s))
    # count number of times each probe is used
    count_check = Counter(probes)
    # make sure that's just once
    assert max(count_check.values()) == 1, "Still some probes in > 1 cluster"


def enter_cluster_as_locus_to_db(c, probe_clusters):
    for cluster, probes in probe_clusters.iteritems():
        for idx in probes:
            c.execute("""UPDATE probes SET locus = ? WHERE id = ?""", (cluster, idx))


def add_new_locus_numbers_to_rest(c, mx):
    c.execute("""SELECT id FROM probes WHERE locus IS NULL""")
    values = c.fetchall()
    for k, idx in enumerate(values):
        new_loc = mx + k
        c.execute("""UPDATE probes SET locus = ? WHERE id = ?""", (new_loc, idx[0]))


def rename_probes_within_each_locus(c):
    c.execute("""SELECT DISTINCT(locus) FROM probes ORDER BY locus""")
    for idx in c.fetchall():
        c.execute("""SELECT id FROM probes WHERE locus = ? ORDER BY id""", idx)
        results = c.fetchall()
        for k, idz in enumerate(results):
            c.execute("""UPDATE probes SET probe = ? WHERE id = ?""", (k + 1, idz[0]))


def main():
    args = get_args()
    # get files
    compare_files = get_compare_files(args)
    # read in data
    taxon_overlaps = get_overlapping_data_by_taxon(compare_files)
    # get all perfect overlaps btw all taxa
    perfect_overlaps = set.intersection(*taxon_overlaps.values())
    # get combinations of 3
    all_overlaps = get_intersection_of_taxa(taxon_overlaps, perfect_overlaps)
    check_probe_counts_in_clusters(all_overlaps)
    probe_clusters = defaultdict(list)
    for k, cluster in enumerate(all_overlaps):
        probe_clusters[k].extend(list(cluster))
    # connect to db
    conn = sqlite3.connect(args.db)
    c = conn.cursor()
    # zero out current locus designations
    c.execute("""UPDATE probes SET locus = Null""")
    enter_cluster_as_locus_to_db(c, probe_clusters)
    add_new_locus_numbers_to_rest(c, max(probe_clusters.keys()))
    # zero out current probe designations
    c.execute("""UPDATE probes SET probe = Null""")
    # rename the probes
    rename_probes_within_each_locus(c)
    conn.commit()
    c.close()
    conn.close()

# locus id 1 and 2


if __name__ == '__main__':
    main()
