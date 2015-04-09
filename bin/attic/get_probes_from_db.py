#!/usr/bin/env python
# encoding: utf-8
"""
File: get_probes_from_db.py
Author: Brant Faircloth

Created by Brant Faircloth on 15 August 2012 14:08 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: 

"""

import argparse
import sqlite3
from phyluce.helpers import is_file, FullPaths

import pdb


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Program description""")
    parser.add_argument(
            "db",
            help="""The database containing probeset tables"""
        )
    parser.add_argument(
            "table",
            type=str,
            help="""The probeset (i.e. column) you want"""
        )
    parser.add_argument(
            "output",
            type=is_file,
            action=FullPaths,
            help="""The file in wich to store the output fasta"""
        )
    parser.add_argument(
            "group",
            type=str,
            choices=['tetrapods','fish'],
            help="""The group to append to the fasta metadata"""
        )
    return parser.parse_args()


def main():
    args = get_args()
    conn = sqlite3.connect(args.db)
    c = conn.cursor()
    # get distinct loci
    query = """SELECT count(distinct(probes.locus))
        FROM probes, probeset
        WHERE probeset.{} = 1
        AND probes.id = probeset.id
        ORDER BY locus, probe""".format(args.table)
    c.execute(query)
    distinct = c.fetchall()[0][0]
    query = """SELECT probes.id, probes.locus, probes.probe, probes.source, probes.sequence
        FROM probes, probeset
        WHERE probeset.{} = 1
        AND probes.id = probeset.id
        ORDER BY locus, probe""".format(args.table)
    #pdb.set_trace()
    c.execute(query)
    probe_data = c.fetchall()
    output = open(args.output, 'w')
    for probe in probe_data:
        idx, locus, probe, source, sequence = probe
        output.write(">uce-{0}_p{1} |source:{2},probes-id:{3},probes-locus:{0},probes-probe:{1},group:{4}\n{5}\n".format(
                locus,
                probe,
                source,
                idx,
                args.group,
                sequence
            ))
    print "{0} probes in set targeting UCE loci".format(len(probe_data))
    print "{0} distinct UCE loci targeted by probes".format(distinct)
    output.close()
    c.close()
    conn.close()

if __name__ == '__main__':
    main()
