#!/usr/bin/env python
# encoding: utf-8
"""
File: get_indels_from_alignments.py
Author: Brant Faircloth

Created by Brant Faircloth on 07 August 2012 21:08 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: 

"""

import os
import sys
import glob
import sqlite3
import argparse
import multiprocessing
from Bio import AlignIO
from Bio.Alphabet import generic_dna
from Bio.Align import MultipleSeqAlignment
from phyluce.helpers import is_dir, FullPaths, get_file_extensions

import pdb


def get_args():
    parser = argparse.ArgumentParser(
            description="""Count indels in alignments, relative to (but excluding) the outgroup"""
        )
    parser.add_argument(
            'input',
            type=is_dir,
            action=FullPaths,
            help="""The directory containing the alignment files"""
        )
    parser.add_argument(
            'outgroup',
            type=str,
            help="""The outgroup taxon"""
        )
    parser.add_argument(
            '--output',
            type=str,
            default='output',
            help="""The output filename (without extension - code will add .sqlite)"""
        )
    parser.add_argument(
            "--input-format",
            dest="input_format",
            choices=['fasta', 'nexus', 'phylip', 'clustal', 'emboss', 'stockholm'],
            default='fasta',
            help="""The input alignment format""",
        )
    parser.add_argument(
            "--trim",
            dest="trim",
            choices=['absolute', 'relative', 'none'],
            default='absolute',
            help="""Count indels only after removing gaps/missing data from alignment ends""",
        )
    parser.add_argument(
            "--cores",
            type=int,
            default=1,
            help="""The number of cores to use.""",
        )
    return parser.parse_args()


def get_files(input_dir, input_format):
    alignments = []
    for ftype in get_file_extensions(input_format):
        alignments.extend(glob.glob(os.path.join(input_dir, "*{}".format(ftype))))
    return alignments


def get_starts_and_stops(seq):
    for position, base in enumerate(seq):
        if base == '-' or base == '?':
            continue
        else:
            start = position
            break
    for position, base in enumerate(seq[::-1]):
        if base == '-' or base == '?':
            continue
        else:
            stop = len(seq) - position
            break
    return start, stop


def get_absolute_starts_and_ends(aln, outgroup):
    """we need to determine actual starts of alignments.  since we're comparing
    to the outgroup, we also need to determine its start and stop, and use those
    values if less than the taxon-specific values"""
    outgroup_start, outgroup_stop = get_starts_and_stops(outgroup)
    trim_positions = {}
    for taxon in aln:
        taxon_start, taxon_stop = get_starts_and_stops(taxon.seq)
        if taxon_start < outgroup_start:
            start = outgroup_start
        else:
            start = taxon_start
        if taxon_stop > outgroup_stop:
            stop = outgroup_stop
        else:
            stop = taxon_stop
        trim_positions[taxon.id] = [start, stop]
    # return the absolute min and max start and stop across all alignments
    minimum = max([v[0] for v in trim_positions.values()])
    maximum = min([v[1] for v in trim_positions.values()])
    trim_positions['minmax'] = [minimum, maximum]
    return trim_positions


def create_indel_database(db):
    """Create the indel database"""
    conn = sqlite3.connect(db)
    c = conn.cursor()
    c.execute("PRAGMA foreign_keys = ON")
    try:
        c.execute('''CREATE TABLE indels (
                idx INTEGER PRIMARY KEY AUTOINCREMENT,
                taxon_name text,
                locus text,
                length float,
                insertions float,
                deletions float,
                substitutions float)'''
            )
    except sqlite3.OperationalError, e:
        if e[0] == 'table indels already exists':
            answer = raw_input("Database already exists.  Overwrite [Y/n]? ")
            if answer == "Y" or "YES":
                os.remove(db)
                conn, c = create_indel_database(db)
            else:
                sys.exit(2)
        else:
            raise sqlite3.OperationalError
            pdb.set_trace()
    return conn, c


def trim_alignments_with_start_and_stop(arguments, iden, new_aln, outgroup, trim_positions):
    if arguments.trim == 'relative':
        trim_new_aln = new_aln[:, trim_positions[iden][0]:trim_positions[iden][1]]
        trim_outgroup = outgroup[trim_positions[iden][0]:trim_positions[iden][1]]
    elif arguments.trim == 'absolute':
        trim_new_aln = new_aln[:, trim_positions['minmax'][0]:trim_positions['minmax'][1]]
        trim_outgroup = outgroup[trim_positions['minmax'][0]:trim_positions['minmax'][1]]
    elif arguments.trim == 'none':
        trim_new_aln = new_aln
        trim_outgroup = outgroup
    return trim_new_aln, trim_outgroup


def worker(work):
    arguments, f = work
    results = {}
    name_map = {}
    try:
        aln = AlignIO.read(f, arguments.input_format)
        # create new alignment to hold everything but the outgroup
        new_aln = MultipleSeqAlignment([], generic_dna)
        count = 0
        for taxon in aln:
            if taxon.id == arguments.outgroup:
                outgroup = taxon
            else:
                results[taxon.id] = {
                        'insertions': [],
                        'deletions': [],
                        'substitutions': [],
                        'length': [],
                        'filename': os.path.splitext(os.path.basename(f))[0]
                    }
                new_aln.append(taxon)
                name_map[taxon.id] = count
                count += 1
        sys.stdout.write('.')
        sys.stdout.flush()
        # exception block catches case where a given alignment has missing data for a taxon
        # which we're avoiding (for now).
        try:
            # get starts and ends of true bases in alignment
            trim_positions = get_absolute_starts_and_ends(new_aln, outgroup)
            # loop through taxa and get outgroup
            bases = set(['A', 'C', 'G', 'T', 'a', 'c', 'g', 't'])
            for taxon in new_aln:
                trim_new_aln, trim_outgroup = trim_alignments_with_start_and_stop(
                            arguments,
                            taxon.id,
                            new_aln,
                            outgroup,
                            trim_positions
                        )
                results[taxon.id]['length'] = trim_new_aln.get_alignment_length()
                for idx in xrange(results[taxon.id]['length']):
                    col = list(trim_new_aln[:, idx])
                    # if outgroup only has insertion it is [ACGT] while set(base) == (['-'])
                    if trim_outgroup[idx] in bases and set(col) == set(['-']):
                        outgroup_insertion = True
                    else:
                        outgroup_insertion = False
                    # if outgroup only has deletion it is '-', while len(set(base)) >= 1 and not (['-'])
                    if trim_outgroup[idx] == '-' and (len(set(col)) == 1 and set(col) != set(['-'])):
                        outgroup_deletion = True
                    else:
                        outgroup_deletion = False
                    if not outgroup_insertion and not outgroup_deletion:
                        # loop through columns and compare to outgroup
                        base = col[name_map[taxon.id]]
                        # insertion
                        if base in bases and trim_outgroup[idx] == '-':
                            results[taxon.id]['insertions'].append(1)
                        # deletion
                        elif base == '-' and trim_outgroup[idx] in bases:
                            results[taxon.id]['deletions'].append(1)
                        # substitution
                        elif base != trim_outgroup[idx]:
                            results[taxon.id]['substitutions'].append(1)
        except UnboundLocalError, e:
            results = {None: os.path.splitext(os.path.basename(f))[0]}
    except ValueError, e:
        if e.message == 'No records found in handle':
            print 'No records found in {0}'.format(os.path.basename(f))
        else:
            raise ValueError('Something is wrong with alignment {0}'.format(os.path.basename(f)))
    return results


def main():
    args = get_args()
    db_name = "{0}-{1}.sqlite".format(args.output, args.trim)
    conn, c = create_indel_database(db_name)
    # iterate through all the files to determine the longest alignment
    work = [(args, f) for f in get_files(args.input, args.input_format)]
    sys.stdout.write("Running")
    if args.cores > 1:
        pool = multiprocessing.Pool(args.cores)
        results = pool.map(worker, work)
    else:
        results = map(worker, work)
    print "\nEntering data to sqlite...."
    for result in results:
        for taxon_name, values in result.iteritems():
            if taxon_name is None:
                print "Did not process alignment {0}. It may have missing data.".format(values)
            else:
                c.execute('''INSERT INTO indels (
                                taxon_name,
                                locus,
                                length,
                                insertions,
                                deletions,
                                substitutions
                            )
                            VALUES (?,?,?,?,?,?)''', (
                                taxon_name,
                                values['filename'],
                                values['length'],
                                sum(values['insertions']),
                                sum(values['deletions']),
                                sum(values['substitutions'])
                            ))
    conn.commit()
    c.close()
    conn.close()
    print "Done."


if __name__ == '__main__':
    main()
