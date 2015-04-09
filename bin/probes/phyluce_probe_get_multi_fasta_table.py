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
from collections import defaultdict
from Bio import SeqIO

from phyluce.helpers import is_dir, is_file, FullPaths
from bx.intervals.intersection import Interval, IntervalTree

import pdb

def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
        description="""Make a table of BED file matches"""
    )
    parser.add_argument(
        "--fastas",
        required=True,
        type=is_dir,
        action=FullPaths,
        help="""A folder of fasta files."""
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
    return parser.parse_args()


def create_match_database(db, organisms, base):
    """Create the UCE-match database"""
    conn = sqlite3.connect(db)
    c = conn.cursor()
    c.execute("PRAGMA foreign_keys = ON")
    try:
        create_string = [org + ' integer DEFAULT 0' for org in organisms]
        query = "CREATE TABLE {0} (locus text primary key, {1})".format(base, ', '.join(create_string))
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
    fastas = glob.glob(os.path.join(args.fastas, "*.fasta"))
    organisms = []
    conserved = defaultdict(list)
    for fasta in fastas:
        taxon_name = os.path.splitext(os.path.basename(fasta))[0]
        organisms.append(taxon_name)
        sys.stdout.write(taxon_name)
        sys.stdout.flush()
        with open(fasta, "rU") as fasta:
            for cnt, seq in enumerate(SeqIO.parse(fasta, "fasta")):
                if cnt % 1000 == 0:
                    sys.stdout.write(".")
                    sys.stdout.flush()
                ids = seq.id.split("|")
                locus = ids[3].split(":")[1]
                conserved[locus].append(taxon_name)
        print ""
    # create database
    print "Creating database"
    conn, c = create_match_database(args.output, organisms, args.base_taxon)
    print "Inserting results"
    for locus, taxa in conserved.iteritems():
        names = ", ".join([i for i in taxa])
        ones = ", ".join(['1'] * len(taxa))
        query = "INSERT INTO {} (locus, {}) values ('{}', {})".format(
            args.base_taxon,
            names,
            locus,
            ones
        )
        c.execute(query)
    conn.commit()
    conn.close()


if __name__ == '__main__':
    main()
