#!/usr/bin/env python
# encoding: utf-8
"""
File: get_probes_from_model_loci_db.py
Author: Brant Faircloth

Created by Brant Faircloth on 08 October 2012 15:10 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: 

"""


import os
import sqlite3
import argparse
from collections import defaultdict

import pdb


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Program description""")
    parser.add_argument(
            "db",
            help="""The database containing the probe matches"""
        )
    parser.add_argument(
            "outfile",
            help="""The output files for the results""",
        )
    return parser.parse_args()


def get_match_data_from_db(cur):
    query = '''SELECT probes.locus, probes.probe, probes.oldlocus, galGal3.*
        FROM probes, probeset, galGal3
        WHERE probes.id = galGal3.id
        AND probes.id = probeset.id
        AND probes5k = 1
        AND source="faircloth"
        '''
    cur.execute(query)
    rows = cur.fetchall()
    return rows


def get_all_locus_coords(rows):
    coords = defaultdict(list)
    chromos = defaultdict(list)
    for row in rows:
        coords[row[2]].append(row[7])
        coords[row[2]].append(row[8])
        chromos[row[2]].append(row[5])
    return chromos, coords


def get_coords(cur, chromos, coords):
    all_data = []
    for locus, coord in coords.iteritems():
        # get coords
        start = min(coord)
        end = max(coord)
        middle = int(round((start + end) / 2, 0))
        # go out 90
        new_start = middle - 91
        new_end = middle + 91
        cur.execute("SELECT locus FROM probes WHERE oldlocus = ? LIMIT 1", (locus,))
        newname = cur.fetchall()[0]
        try:
            assert len(set(chromos[locus])) == 1
            all_data.append({'name':newname[0], 'chromo':chromos[locus][0], 'start':new_start, 'end':new_end, 'oldname':locus})
        except:
            print "Skipping uce-{}, oldname = {}".format(newname[0], locus)
    return all_data



def main():
    args = get_args()
    conn = sqlite3.connect(args.db)
    cur = conn.cursor()
    rows = get_match_data_from_db(cur)
    chromos, coords = get_all_locus_coords(rows)
    all_data = get_coords(cur, chromos, coords)
    outf = open(args.outfile, 'w')
    outf.write("locus,chromo,new-start,new-end,old-name,map-position\n")
    #pdb.set_trace()
    for row in all_data:
        outf.write("uce-{0[name]},{0[chromo]},{0[start]},{0[end]},{0[oldname]},{0[chromo]}:{0[start]}-{0[end]}\n".format(
            row
            ))
    outf.close()
    outf = open(args.outfile + ".bed", 'w')
    outf.write('''track name=uces description="UCEs with new coordinates" useScore=0\n''')
    for row in all_data:
        outf.write("{0[chromo]}\t{0[start]}\t{0[end]}\tuce-{0[name]}\n".format(
            row
            ))


if __name__ == '__main__':
    main()