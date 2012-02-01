#!/usr/bin/env python
# encoding: utf-8

"""
get_sureselect_probes_or_cons.py

Created by Brant Faircloth on 31 May 2011.
Copyright 2011 Brant C. Faircloth. All rights reserved.

Given an input database (probes.sqlite) get the sequence of either the 
probes or the UCEs that the probes target.

Also, given the --genome option, generate a fake genome composed of
the probes that target each, respective UCE.  In many cases, this is
only one probe.  In others, this is multiple probes, so we have to join
them together to get our representation of the UCE (which is actually the
UCE ± a little but of flank).  After doing this, write these out, placing
each uce region 1 kb from the next, and record a metadata file giving the
positions of everything, so we know what's what.

"""

import os
import sys
import sqlite3
import argparse
import operator
import textwrap
from collections import defaultdict


import pdb

def get_args():
    parser = argparse.ArgumentParser(description='Get UCE or probe sequence from a database')
    parser.add_argument('database', help = 'The database to query')
    parser.add_argument('--output', nargs='?', type=argparse.FileType('w'), default=sys.stdout)
    parser.add_argument('--genome', help = 'Generate fake genome', action='store_true')
    parser.add_argument('--all', help = 'Return all conservation probes', action='store_false')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--conserved', action='store_true')
    group.add_argument('--probes', action='store_false')
    return parser.parse_args()

def get_conserved_sequences(cur):
    """docstring for get_conserved_sequences"""
    cur.execute("SELECT sureselect_probe_counts.'sureselect.seq', \
        sureselect_probe_counts.cnt, \
        sureselect_probe_counts.data_source, \
        cons.cons \
        FROM sureselect_probe_counts, cons \
        WHERE sureselect_probe_counts.'sureselect.id' = cons.id \
        AND data_source = 'conservation'")
    return cur.fetchall()

def get_probe_sequences(cur, selected = True):
    """docstring for get_probe_sequences"""
    if selected:
        cur.execute("SELECT sureselect_probe_counts.'sureselect.seq', \
            sureselect_probe_counts.cnt, \
            sureselect_probe_counts.data_source, \
            sureselect.probe_name, \
            sureselect.probe_sequence \
            FROM sureselect_probe_counts, sureselect \
            WHERE sureselect_probe_counts.'sureselect.id' = sureselect.id \
            AND sureselect.data_source = 'conservation' \
            AND sureselect.selected = 1")
    else:
        cur.execute("SELECT sureselect_probe_counts.'sureselect.seq', \
            sureselect_probe_counts.cnt, \
            sureselect_probe_counts.data_source, \
            sureselect.probe_name, \
            sureselect.probe_sequence \
            FROM sureselect_probe_counts, sureselect \
            WHERE sureselect_probe_counts.'sureselect.id' = sureselect.id \
            AND sureselect.data_source = 'conservation'")
    return cur.fetchall()

def main():
    args = get_args()
    conn = sqlite3.connect(args.database)
    cur = conn.cursor()
    if args.conserved:
        seq = get_conserved_sequences(cur)
    else:
        seq = get_probe_sequences(cur, args.all)
    if not args.genome:
        for s in seq:
            header = ">{0}|{1}|probes:{2}".format(s[0].lower(), s[2], s[1])
            txt = textwrap.wrap(s[-1], 78)
            args.output.write("{0}\n{1}\n".format(header, '\n'.join(txt)))
    elif not args.conserved and args.genome:
        genome = []
        genome.append('N' * 1000)
        # really need to assemble probe sequences and then space those out in a
        # pseudo-genome
        cons = defaultdict(list)
        for s in seq:
            ss = s[3].split('_')
            name = s[0]
            if name.startswith('chrE22') or name.startswith('chrUn'):
                cons[name].append([ss[-1], '-'.join([ss[3],ss[4]]), s[-1]])
            else:
                #pdb.set_trace()
                cons[name].append([ss[-1], ss[2], s[-1]])                
        count = 1
        metadata = {}
        for k,v in sorted(cons.iteritems(), key=operator.itemgetter(0)):
            #pdb.set_trace()
            if len(v) == 1:
                target_map = v[0][1]
                target = v[0][-1]
            else:
                bits = []
                for pos, probe in enumerate(v):
                    #pdb.set_trace()
                    assert pos == int(probe[0]), "positions don't match"
                    if pos == 0:
                        bracket = [probe[1]]
                    if pos + 1 < len(v):
                        start = int(probe[1].split(':')[1].split('-')[0])
                        next_start = int(v[pos + 1][1].split(':')[1].split('-')[0])
                        size = next_start - start
                        bits.append(probe[2][:size])
                    else:
                        bits.append(probe[2])
                        bracket.append(probe[1])
                target_map = '-'.join([bracket[0].split('-')[0], bracket[1].split('-')[1]])
                target = ''.join(bits)
            try:
                start, stop = target_map.split(':')[1].split('-')
            except:
                pdb.set_trace()
            length = int(stop) - int(start) + 1
            metadata[count * 1000] = [k, length, target_map]
            count += 1
            genome.append(target + 'N' * (1000 - len(target)))
        synth_genome = ''.join(genome)
        args.output.write('>UCEProbes Synthetic Genome\n{0}'.format(synth_genome))
        #pdb.set_trace()
        os.system('fold -w 60 < {0} > {0}.folded'.format(args.output.name))
        #pdb.set_trace()
        mk = metadata.keys()
        mk.sort()
        mout = open('{0}.metadata.csv'.format(args.output.name), 'w')
        for k in mk:
            mout.write("{0},{1},{2},{3}\n".format(metadata[k][0], k, metadata[k][1], metadata[k][2]))
        mout.close()
    

if __name__ == '__main__':
    main()
