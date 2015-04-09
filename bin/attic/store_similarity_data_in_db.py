#!/usr/bin/env python
# encoding: utf-8
"""
File: store_similarity_data_in_db.py
Author: Brant Faircloth

Created by Brant Faircloth on 07 April 2012 11:04 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: 

"""

import os
import re
import sys
import sqlite3
import argparse
import ConfigParser
from collections import defaultdict

from phyluce.helpers import FullPaths, is_file, get_uce_names_from_probes

#import pdb


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Given a config file of similarity results, enter those to database""")
    parser.add_argument(
            "conf",
            action=FullPaths,
            type=is_file,
            help="""The path to the config file containing directory mappings"""
        )
    parser.add_argument(
            "probes",
            action=FullPaths,
            type=is_file,
            help='Fasta file of the probe sequence searched',
        )
    parser.add_argument(
            "output",
            action=FullPaths,
            help="""The name of the output database""",
        )
    return parser.parse_args()


def create_similarity_database(db, organisms, uces):
    """docstring for create_probe_database"""
    conn = sqlite3.connect(db)
    c = conn.cursor()
    c.execute("PRAGMA foreign_keys = ON")
    try:
        tables = ['diverge_gaps', 'diverge_no_gaps', 'distance', 'length']
        create_string = [org + ' float' for org in organisms]
        for table in tables:
            query = "CREATE TABLE {0} (uce text primary key, {1})".format(
                    table,
                    ','.join(create_string)
                )
            c.execute(query)
        for table in tables:
            for uce in uces:
                query = "INSERT INTO {0}(uce) values ('{1}')".format(table, uce)
                c.execute(query)
    except sqlite3.OperationalError, e:
        if e[0] == 'table diverge_gaps already exists':
            answer = raw_input("Database already exists.  Overwrite [Y/n]? ")
            if answer == "Y" or "YES":
                os.remove(db)
                conn, c = create_similarity_database(db, organisms, uces)
            else:
                sys.exit(2)
        else:
            raise sqlite3.OperationalError
    return conn, c


def get_organism_name_from_conf(filemap, conf):
    return [element[0] for element in filemap]


def get_data_from_files(directory, files):
    pdd, sdd = defaultdict(list), defaultdict(list)
    for line in open(os.path.join(directory, files[0]), 'rU'):
        if not line.startswith('locus'):
            ls = line.strip().split(',')
            pdd[ls[0]] = ls[1:]
    for line in open(os.path.join(directory, files[1]), 'rU'):
        if not line.startswith('locus'):
            ls = line.strip().split(',')
            sdd[ls[0]] = ls[1:]
    return pdd, sdd


def main():
    args = get_args()
    conf = ConfigParser.ConfigParser()
    conf.read(args.conf)
    filemap = conf.items('taxa')
    regex = re.compile("s_[0-9]+$")
    organisms = get_organism_name_from_conf(filemap, conf)
    #pdb.set_trace()
    # get list of all loci to create tables
    uces = get_uce_names_from_probes(args.probes, regex, 's')
    # create db and return conn and cur
    conn, cur = create_similarity_database(args.output, organisms, uces)
    # iterate over output files to grab values and store in db
    files = ['pairwise_distance.csv', 'sequence_divergence.csv']
    for values in filemap:
        organism, directory = values
        print "Processing {0}...".format(organism)
        sd, pd = get_data_from_files(directory, files)
        for locus, values in sd.iteritems():
            locus = locus.strip('.nex')
            query = "UPDATE distance SET {0}={1} WHERE uce = '{2}'".format(organism, values[0], locus)
            cur.execute(query)
        for locus, values in pd.iteritems():
            locus = locus.strip('.nex')
            query = "UPDATE length SET {0}={1} WHERE uce = '{2}'".format(organism, values[0], locus)
            cur.execute(query)
            query = "UPDATE diverge_gaps SET {0}={1} WHERE uce = '{2}'".format(organism, values[1], locus)
            cur.execute(query)
            query = "UPDATE diverge_no_gaps SET {0}={1} WHERE uce = '{2}'".format(organism, values[2], locus)
            cur.execute(query)
    conn.commit()
    cur.close()
    conn.close()


if __name__ == '__main__':
    main()
