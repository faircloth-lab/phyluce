#!/usr/bin/env python
# encoding: utf-8

"""
get_probes_matches_from_lastz.py

Created by Brant Faircloth on 23 June 2011.
Copyright 2011 Brant C. Faircloth. All rights reserved.
"""


import os
import sys
import glob
import sqlite3
import argparse
import textwrap
import bx.seq.twobit
from operator import itemgetter
from collections import defaultdict
from seqtools.sequence import fasta, transform
from phyluce import lastz
from phyluce.helpers import get_name, get_dupes, get_matches, run_checks

import pdb

def is_dir(dirname):
    if not os.path.isdir:
        msg = "{0} is not a directory".format(dirname)
        raise argparse.ArgumentTypeError(msg)
    else:
        return dirname

def get_args():
    parser = argparse.ArgumentParser(description='Match UCE probes to assembled contigs and store the data')
    parser.add_argument('lastz', help='A directory of lastz output', type=is_dir)
    parser.add_argument('query', help='The fasta or 2bit file containing the probe sequence')
    parser.add_argument('db', help='The name of the output db')
    parser.add_argument('--name-components', dest = 'components', help = 'number of parts in the name', default = 2, type = int)
    parser.add_argument('--splitchar', help = 'The name character on which to split', default = "_", type = str)
    parser.add_argument('--dupefile', help='The path to a lastz file of lastz-against-self results')
    parser.add_argument('--verbose', action='store_true', default = False)
    return parser.parse_args()

def create_match_database(db, organisms, uces):
    """docstring for create_probe_database"""
    conn = sqlite3.connect(db)
    c = conn.cursor()
    c.execute("PRAGMA foreign_keys = ON")
    try:
        #pdb.set_trace()
        create_string = [org + ' text' for org in organisms]
        query = "CREATE TABLE matches (uce text primary key, {0})".format(','.join(create_string))
        c.execute(query)
        for uce in uces:
            c.execute("INSERT INTO matches(uce) values (?)", (uce,))
    except sqlite3.OperationalError, e:
        if e[0] == 'table matches already exists':
            answer = raw_input("Database already exists.  Overwrite [Y/n]? ")
            if answer == "Y" or "YES":
                os.remove(db)
                conn, c = create_match_database(db, organisms, uces)
            else:
                sys.exit(2)
        else:
            raise sqlite3.OperationalError
    return conn, c

def store_lastz_results_in_db(c, critter, region):
    insert_string = "UPDATE matches SET {0} = 1 WHERE uce = '{1}'".format(critter, region)
    c.execute(insert_string)

def main():
    args = get_args()
    uces = set([get_name(read.identifier, "|", 1) for read in fasta.FastaReader(args.query)])
    files = glob.glob(os.path.join(args.lastz, '*.lastz'))
    # this prob. needs to be more robust
    organisms = [os.path.splitext(os.path.basename(f).split('-')[-1])[0].replace('-',"_") for f in files]
    conn, c = create_match_database(args.db, organisms, uces)
    if args.dupefile:
        dupes = get_dupes(args.dupefile)
    else:
        dupes = None
    #pdb.set_trace()
    for f in files:
        critter = os.path.splitext(os.path.basename(f).split('-')[-1])[0]
        matches, probes = get_matches(f, args.splitchar, args.components)
        count = 0
        for k,v in matches.iteritems():
            skip = False
            if len(v) > 1:
                if run_checks(k, v, probes, args.verbose):
                    # sort by match position
                    v_sort = sorted(v, key = itemgetter(2))
                    start, end = v_sort[0][2], v_sort[-1][3]
                    diff = end - start
                    # ensure our range is less than N(probes) * probe_length - this
                    # still gives us a little wiggle room because probes are ~ 2X tiled
                    if diff > (probes[k] * 120):
                        skip = True
                        if args.verbose:
                            print "range longer than expected"
                else:
                    skip = True
            elif args.dupefile and k in dupes:
                skip = True
                if args.verbose:print "{0} is in dupefile".format(k)
            else:
                pass
            if not skip:
                store_lastz_results_in_db(c, critter, k)
                count += 1
        print "Entered {} matches for {}".format(count, critter)
    conn.commit()
    c.close()
    conn.close()

if __name__ == '__main__':
    main()
