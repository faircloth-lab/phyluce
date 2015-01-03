#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2014 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 11 April 2014 11:41 PDT (-0700)
"""


import os
import sys
import glob
import sqlite3
import argparse
import ConfigParser
from collections import defaultdict

from phyluce.helpers import is_dir, is_file, FullPaths
from bx.intervals.intersection import Interval, IntervalTree

import pdb

def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
        description="""Make a table of BED file matches"""
    )
    parser.add_argument(
        "--conf",
        required=True,
        type=is_file,
        action=FullPaths,
        help="""A config file mapping names to BED files."""
    )
    parser.add_argument(
        "--output",
        required=True,
        help="""A SQLite database to create during integration."""
    )
    parser.add_argument(
        "--base-taxon",
        required=True,
        type=str,
        help="""The base taxon to use."""
    )
    parser.add_argument(
        "--filter",
        type=str,
        help="""A file-type filter to apply to the BED directory of files""",
    )
    return parser.parse_args()


def create_match_database(db, organisms, base):
    """Create the UCE-match database"""
    conn = sqlite3.connect(db)
    c = conn.cursor()
    c.execute("PRAGMA foreign_keys = ON")
    try:
        create_string = [org + ' integer DEFAULT 0' for org in organisms]
        query = "CREATE TABLE {0} (uce integer primary key autoincrement, chromo text, start integer, stop integer, {1})".format(base, ', '.join(create_string))
        c.execute(query)
    except sqlite3.OperationalError, e:
        #log.critical("Database already exists")
        if e[0].endswith("already exists"):
            answer = raw_input("Database already exists.  Overwrite [Y/n]? ")
            if answer == "Y" or "YES":
                os.remove(db)
                conn, c = create_match_database(db, organisms, base)
            else:
                sys.exit(2)
        else:
            #log.critical("Cannot create database")
            raise sqlite3.OperationalError("Cannot create database")
    return conn, c


def main():
    args = get_args()
    conf = ConfigParser.ConfigParser()
    conf.read(args.conf)
    organisms = []
    beds = conf.items("beds")
    conserved = defaultdict(IntervalTree)
    for taxon_name, bedfile in beds:
        organisms.append(taxon_name)
        sys.stdout.write(taxon_name)
        sys.stdout.flush()
        with open(bedfile, "rU") as bed:
            for cnt, line in enumerate(bed):
                if cnt % 1000 == 0:
                    sys.stdout.write(".")
                    sys.stdout.flush()
                chromo, start, stop = [int(i) if i.isdigit() else i for i in line.strip().split('\t')]
                if chromo not in conserved.keys():
                    conserved[chromo].insert_interval(Interval(start, stop, set([taxon_name])))
                elif chromo in conserved.keys():
                    # check for intersection
                    overlaps = conserved[chromo].find(start, stop)
                    if overlaps:
                        for o in overlaps:
                                o.value.add(taxon_name)
                    else:
                        conserved[chromo].insert_interval(Interval(start, stop, set([taxon_name])))
        print ""
    # create database
    print "Creating database"
    conn, c = create_match_database(args.output, organisms, args.base_taxon)
    print "Inserting results"
    for chromo, ints in conserved.iteritems():
        for each in ints.after(0,1000000000,1000000000):
            names = ", ".join([i for i in list(each.value)])
            ones = ", ".join(['1'] * len(each.value))
            query = "INSERT INTO {}(chromo, start, stop, {}) values ('{}', {}, {}, {})".format(
                args.base_taxon,
                names,
                chromo,
                each.start,
                each.end,
                ones
            )
            c.execute(query)
    conn.commit()
    conn.close()


if __name__ == '__main__':
    main()
