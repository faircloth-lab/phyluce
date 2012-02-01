#!/usr/bin/env python
# encoding: utf-8

"""
get_match_counts.py

Created by Brant Faircloth on 03 June 2011.
Copyright 2011 Brant C. Faircloth. All rights reserved.

Used to determine, for a group specified within a config file, the number of
UCE loci returned for all members of that group.  Also returns optimum group
sizes based on complete enumeration and sampling approaces.


"""

import os
import sys
import random
import sqlite3
import operator
import argparse
import itertools
import ConfigParser
import multiprocessing
from collections import Counter
from collections import defaultdict

import pdb

def get_args():
    parser = argparse.ArgumentParser(description='Match UCE probes to assembled contigs and store the data')
    parser.add_argument('db', help='The database holding probe matches')
    parser.add_argument('config', help='The config file containing match information')
    parser.add_argument('group', help='The config group whose results you want', type=str)
    parser.add_argument('--output', nargs='?', type=argparse.FileType('w'), default=sys.stdout)
    parser.add_argument('--probes', help='Return data on probles only', action="store_true")
    parser.add_argument('--optimize', help='Return optimum groups of probes by enumeration (default) or sampling.', action="store_true")
    parser.add_argument('--random', help='Optimize by sampling', action="store_true")
    parser.add_argument('--samples', help='The number of samples to take', type = int, default = 10)
    parser.add_argument('--sample-size', dest = 'sample_size', help='The group size of samples', type = int, default = 10)
    parser.add_argument('--extend', dest = 'extend', help='The match database to add as an extension')
    parser.add_argument('--notstrict', dest = 'notstrict', help='Do not do strict matching', action="store_true")
    parser.add_argument('--silent', dest = 'silent', help='Don\'t print probe names', action="store_true")
    parser.add_argument('--keep-counts', dest = 'keep_counts', action = "store_true")
    return parser.parse_args()

def get_name(header, splitchar = "_", items = 2):
    return "_".join(header.split(splitchar)[:items]).lstrip('>')

def get_uce_names(c):
    c.execute("SELECT uce FROM matches")
    return set([uce[0] for uce in c.fetchall()])

def get_names_from_config(config, group):
    try:
        return [i[0] for i in config.items(group)]
    except ConfigParser.NoSectionError:
        return None

def remove_duplicates_from(c, organism, matches):
    if not organism.endswith('*'):
        query = "SELECT uce, {0} FROM match_map WHERE uce in ({1})".format(organism, ','.join(["'{0}'".format(i) for i in matches]))
    else:
        query = "SELECT uce, {0} FROM extended.match_map WHERE uce in ({1})".format(organism.rstrip('*'), ','.join(["'{0}'".format(i) for i in matches]))
    c.execute(query)
    data = c.fetchall()
    m = defaultdict(list)
    for d in data:
        node = d[1].split('(')[0]
        m[node].append(d[0])
    single_matches = [v[0] for k,v in m.iteritems() if len(v) <= 1]
    #if 'monodelphis_domestica' in organism:
    #    pdb.set_trace()
    return single_matches

def get_all_matches_by_organism(c, organisms):
    organismal_matches = {}
    for organism in organisms:
        if not organism.endswith('*'):
            c.execute("SELECT uce FROM matches WHERE {0} = 1".format(organism))
        else:
            c.execute("SELECT uce FROM extended.matches WHERE {0} = 1".format(organism.rstrip('*')))
        matches = set([uce[0] for uce in c.fetchall()])
        # we've removed dupe UCE matches, but we need to remove
        # dupe node matches (i,e. we've removed dupe target matches
        # and we also need to remove dupe query matches - the pop up as
        # data change, so it's a constant battle)
        matches = remove_duplicates_from(c, organism, matches)
        organismal_matches[organism.rstrip('*')] = matches
    return organismal_matches

def return_group_matches(organismal_matches, organisms, uces, fast = True):
    losses = {}
    if fast:
        setlist = [organismal_matches[organism.rstrip('*')] for organism in organisms]
        # this is faster than loop (by 50%)
        uces = uces.intersection(*setlist)
    else:
        for organism in organisms:
            old_uce_len = len(uces)
            uces = uces.intersection(organismal_matches[organism.rstrip('*')])
            losses[organism.rstrip('*')] = old_uce_len - len(uces)
    return uces, losses

def return_notstrict_matches(organismal_matches, organisms, uces):
    setlist = [organismal_matches[organism.rstrip('*')] for organism in organisms]
    return uces.union(*setlist), None

def optimize_group_match_runner(combos):
    organismal_matches, organisms, size, uces = combos
    mx = None
    for group in itertools.combinations(organisms, size):
        group_uces, losses = return_group_matches(organismal_matches, group, uces)
        if len(group_uces) > mx:
            best = group
            mx = len(group_uces)
            best_uces = group_uces
    sys.stdout.write(".")
    sys.stdout.flush()
    return [best, mx, best_uces]

def optimize_group_matches(c, organisms, uces, random):
    combos = []
    organismal_matches = get_all_matches_by_organism(c, organisms)
    if not random:
        for size in xrange(1, len(organisms) + 1):
            combos.append([organismal_matches, organisms, size, uces])
    else:
        combos.append([organismal_matches, organisms, len(organisms) - 1, uces])
    if not random:
        sys.stdout.write("Processing")
        sys.stdout.flush()
        pool = multiprocessing.Pool(6)
        results = pool.map(optimize_group_match_runner, combos)
    else:
        results = map(optimize_group_match_runner, combos)
    return results
    
def main():
    args = get_args()
    config = ConfigParser.RawConfigParser(allow_no_value=True)
    config.read(args.config)
    conn = sqlite3.connect(args.db)
    c = conn.cursor()
    if args.extend:
        query = "ATTACH DATABASE '{0}' AS extended".format(args.extend)
        c.execute(query)
    organisms = get_names_from_config(config, args.group)
    excludes = get_names_from_config(config, 'Excludes')
    if excludes:
        organisms = [org for org in organisms if org not in excludes]
    uces = get_uce_names(c)
    if args.keep_counts:
        all_counts = []
    if args.optimize:
        if not args.random:
            # this is the complete enumeration approach - slow
            # for large groups.  Also used mostly to show
            # how groups increase.  Not as powerful as sampling.
            results = optimize_group_matches(c, organisms, uces, False)
            print ""
            for r in results:
                print "{0}\t{1}".format(','.join(r[0]), r[1])
            best_group = None
        else:
            mx = None
            sys.stdout.write("Sampling ")
            sys.stdout.flush()
            missing = []
            for i in xrange(args.samples):
                sys.stdout.write(".")
                sys.stdout.flush()
                # create groups of sample size + 1 so we can look at all
                # combinations of desired sample size (there is only 1
                # combination of all individuals)
                orgs = random.sample(organisms, args.sample_size + 1)
                group, group_size, group_uces = optimize_group_matches(c, orgs, uces, True)[0]
                #sys.stdout.write("[{0}] ".format(group_size))
                #sys.stdout.flush()
                missing.extend([i for i in set(organisms).difference(set(group))])
                if group_size > mx:
                    mx = group_size
                    best_group = group
                    best_uces = group_uces
                if args.keep_counts:
                    all_counts.append((args.sample_size, group_size))
            if not args.keep_counts:
                print "\nmax UCE = {0}".format(mx)
                print "group size = {0}".format(len(best_group))
                print "best group\n\t{0}\n".format(sorted(best_group))
                print "Times not in best group per iteration\n\t{0}\n".format(Counter(missing))
            else:
                args.output.write('\n'.join(["{},{}".format(str(i), str(j)) for
                    i,j in all_counts]))
    else:
        organismal_matches = get_all_matches_by_organism(c, organisms)
        #pdb.set_trace()
        if not args.notstrict:
            best_uces, losses = return_group_matches(organismal_matches, organisms, uces, fast = False)
        else:
            best_uces, losses = return_notstrict_matches(organismal_matches, organisms, uces)
        best_group = organisms
        print "Shared UCEs: {0}\n".format(len(best_uces))
        if losses:
            sorted_losses = sorted(losses.iteritems(), key=operator.itemgetter(1))
            sorted_losses.reverse()
            for loss in sorted_losses:
                print "{0}:{1}".format(loss[0], loss[1])
    if args.output and best_group and not args.silent and not args.keep_counts:
        args.output.write("[Organisms]\n{0}\n[Loci]\n{1}\n".format(
                    '\n'.join(sorted(best_group)), 
                    '\n'.join(sorted(best_uces))))
    if args.probes:
        for uce in get_uce_names(c):
            query = "SELECT {0} FROM matches WHERE uce = '{1}'".format(', '.join(organisms), uce)
            c.execute(query)
            row = c.fetchall() 
            presence = [val for val in row[0][1:]]
            #pdb.set_trace()
            missing = (1 - float(sum(presence))/len(presence)) * 100
            print "{0},{1},{2},{3:.2f}%".format(uce, presence.count(1), presence.count(0), missing)
    args.output.close()
    #pdb.set_trace()

if __name__ == '__main__':
    main()
