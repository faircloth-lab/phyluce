#!/usr/bin/env python
# encoding: utf-8

"""
match_contigs_to_probes.py

Created by Brant Faircloth on 02 June 2011.
Copyright 2011 Brant C. Faircloth. All rights reserved.
"""

import os
import sys
import glob
import copy
import sqlite3
import argparse
from phyluce import lastz
from collections import defaultdict
from seqtools.sequence import fasta

import pdb

def is_dir(dirname):
    if not os.path.isdir:
        msg = "{0} is not a directory".format(dirname)
        raise argparse.ArgumentTypeError(msg)
    else:
        return dirname

def get_args():
    parser = argparse.ArgumentParser(description='Match UCE probes to assembled contigs and store the data')
    parser.add_argument('contigs', help='The directory containing the contigs to match against probes', type=is_dir)
    parser.add_argument('query', help='The query fasta or 2bit file')
    parser.add_argument('output', help='The directory in which to store the lastz alignments', type=is_dir)
    parser.add_argument('--coverage', default = 80, type = int)
    parser.add_argument('--identity', default = 80, type = int)
    parser.add_argument('--dupefile', help='Path to self-to-self lastz results')
    
    return parser.parse_args()

def get_summary_stats_from_parsed_lastz():
    """use summary values from spreadsheets as stopgap,
    but soon, paste count, trim, assembly, etc data into
    a config file so we can easily parse out whatever values
    we need for downstream operations (e.g. count of contigs,
    etc.)
    
    singleton matches
    doublet matches
    triplet matches
    average percent identity
    average percent coverage
    percent contigs matching
    
    """
    pass

def create_probe_database(db, organisms, uces):
    """docstring for create_probe_database"""
    conn = sqlite3.connect(db)
    c = conn.cursor()
    c.execute("PRAGMA foreign_keys = ON")
    try:
        create_string = [org + ' text' for org in organisms]
        query = "CREATE TABLE matches (uce text primary key, {0})".format(','.join(create_string))
        c.execute(query)
        query = "CREATE TABLE match_map (uce text primary key, {0})".format(','.join(create_string))
        c.execute(query)
        for uce in uces:
            c.execute("INSERT INTO matches(uce) values (?)", (uce,))
            c.execute("INSERT INTO match_map(uce) values (?)", (uce,))
    except sqlite3.OperationalError, e:
        if e[0] == 'table matches already exists':
            answer = raw_input("Database already exists.  Overwrite [Y/n]? ")
            if answer == "Y" or "YES":
                os.remove(db)
                conn, c = create_probe_database(db, organisms, uces)
            else:
                sys.exit(2)
        else:
            raise sqlite3.OperationalError 
            pdb.set_trace()
    return conn, c

def store_lastz_results_in_db(c, matches, orientation, critter):
    count = 0
    for key, match in matches.iteritems():
        # we need to check probes for dupe matches to multiple nodes
        if len(match) == 1:
            insert_string = "UPDATE matches SET {0} = 1 WHERE uce = '{1}'".format(critter, list(match)[0])
            c.execute(insert_string)
            #pdb.set_trace()
            orient_key = "{0}({1})".format(key, list(orientation[list(match)[0]])[0])
            insert_string = "UPDATE match_map SET {0} = '{1}' WHERE uce = '{2}'".format(critter, orient_key, list(match)[0])
            c.execute(insert_string)
        else:
            count += 1
    return count

def get_name(header, splitchar = "_", items = 2):
    return "_".join(header.split(splitchar)[:items]).lstrip('>').lower()

def get_dupes(lastz_file):
    matches = defaultdict(list)
    dupes = set()
    for lz in lastz.Reader(lastz_file):
        target_name = get_name(lz.name1, "|", 1)
        query_name = get_name(lz.name2, "|", 1)
        matches[target_name].append(query_name)
    # see if one probe matches any other probes
    # other than the children of the locus
    for k,v in matches.iteritems():
        # if the probe doesn't match itself, we have
        # problems
        if len(v) > 1:
            for i in v:
                if i != k:
                    dupes.add(k)
                    dupes.add(i)
        elif k != v[0]:
            dupes.add(k)
    return dupes

def main():
    args = get_args()
    files = glob.glob(os.path.join(args.contigs, '*.fa*'))
    organisms = [os.path.basename(f).split('.')[0].replace('-',"_") for f in files]
    # make a set to get uniques
    uces = set([get_name(read.identifier, "|", 1) for read in fasta.FastaReader(args.query)])
    conn, c = create_probe_database(os.path.join(args.output, 'probe.matches.sqlite'), organisms, uces)
    #pdb.set_trace()
    print "Processing:"
    if args.dupefile:
        print "\t Getting dupes"
        dupes = get_dupes(args.dupefile)
    for contig in files:
        critter = os.path.basename(contig).split('.')[0].replace('-',"_")
        output = os.path.join(args.output, 
                    os.path.splitext(os.path.basename(contig))[0] + '.lastz')
        contigs = sum([1 for line in open(contig, 'rU').readlines() if line.startswith('>')])
        # align the probes to the contigs
        alignment = lastz.Align(contig, args.query, args.coverage, \
                            args.identity, output)
        lzstdout, lztstderr = alignment.run()
        # parse the lastz results of the alignment
        matches = defaultdict(set)
        orientation = defaultdict(set)
        revmatch = defaultdict(set)
        probe_dupes = set()
        if not lztstderr:
            for lz in lastz.Reader(output):
                # get strandedness of match
                contig_name = get_name(lz.name1) #+ "({0})".format(lz.strand2)
                uce_name = get_name(lz.name2, "|", 1)
                if args.dupefile and uce_name in dupes:
                    probe_dupes.add(uce_name)
                else:
                    matches[contig_name].add(uce_name)
                    orientation[uce_name].add(lz.strand2)
                    revmatch[uce_name].add(contig_name)
        # we need to check nodes for dupe matches to the same probes
        test = []
        node_dupe_matches = set([i for uce,node in revmatch.iteritems() if len(node) > 1 for i in list(node)])
        # pull out our dupe and/or dubious nodes
        oldm = copy.deepcopy(matches)
        for k in matches.keys():
            if k in node_dupe_matches:
                matches.pop(k)
        probe_dupe_matches = store_lastz_results_in_db(c, matches, orientation, critter)
        #pdb.set_trace()
        conn.commit()
        unique_matches = sum([1 for node, uce in matches.iteritems() if len(uce) <= 1])
        out = "\t {0}: {1} ({2:.2f}%) uniques of {3} contigs, {4} dupe probes, "+ \
            "{5} dupe probe matches, {6} dupe node matches"
        print out.format(critter, unique_matches, float(unique_matches)/contigs * 100, contigs, 
                    len(probe_dupes), probe_dupe_matches, len(node_dupe_matches))
        #pdb.set_trace()
            

if __name__ == '__main__':
    main()
