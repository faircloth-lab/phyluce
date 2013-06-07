#!/usr/bin/env python
# encoding: utf-8

"""
match_contigs_to_probes.py

Created by Brant Faircloth on 02 June 2011.
Copyright 2011 Brant C. Faircloth. All rights reserved.
"""

import re
import os
import sys
import glob
import copy
import sqlite3
import argparse
from phyluce import lastz
from phyluce.helpers import is_dir, is_file
from collections import defaultdict
from seqtools.sequence import fasta

import pdb



def get_args():
    parser = argparse.ArgumentParser(description='Match UCE probes to assembled contigs and store the data')
    parser.add_argument(
            'contigs',
            type=is_dir,
            help="The directory containing the contigs to match against probes"
        )
    parser.add_argument('query',
            type=is_file,
            help="The query fasta or 2bit file"
        )
    parser.add_argument(
            'output',
            type=is_dir,
            help="The directory in which to store the lastz alignments"
        )
    parser.add_argument(
            '--coverage',
            default=80,
            type=int)
    parser.add_argument(
            '--identity',
            default=80,
            type=int)
    parser.add_argument(
            '--dupefile',
            help="Path to self-to-self lastz results"
        )
    parser.add_argument(
            "--regex",
            type=str,
            default=None,
            help="""A regular expression to apply to the probe sequences for replacement""",
        )
    parser.add_argument(
            "--repl",
            type=str,
            default=None,
            help="""The replacement text for matches to the regular expression in --regex""",
        )
    args = parser.parse_args()
    if args.regex is not None and args.repl is None:
        sys.exit("If you are replacing text with a regular expression you must pass args.repl value")
    elif args.repl is not None and args.regex is None:
        sys.exit("If you are replacing text with a regular expression you must pass args.regex value")
    return args


def create_probe_database(db, organisms, uces):
    """Create the UCE-match database"""
    conn = sqlite3.connect(db)
    c = conn.cursor()
    c.execute("PRAGMA foreign_keys = ON")
    try:
        create_string = [org + ' text' for org in organisms]
        query = "CREATE TABLE matches (uce text primary key, {0})".format(','.join(create_string))
        c.execute(query)
        query = "CREATE TABLE match_map (uce text primary key, {0})".format(','.join(create_string))
        c.execute(query)
        # convert uces to list of tuples for executemany
        all_uces = [(uce,) for uce in uces]
        c.executemany("INSERT INTO matches(uce) values (?)", all_uces)
        c.executemany("INSERT INTO match_map(uce) values (?)", all_uces)
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
    """enter matched loci in database"""
    for key, match in matches.iteritems():
        # We should have dropped all duplicates at this point
        assert len(match) == 1, "More than one match"
        item = list(match)[0]
        insert_string = "UPDATE matches SET {0} = 1 WHERE uce = '{1}'".format(critter, item)
        c.execute(insert_string)
        #pdb.set_trace()
        orient_key = "{0}({1})".format(key, list(orientation[item])[0])
        insert_string = "UPDATE match_map SET {0} = '{1}' WHERE uce = '{2}'".format(critter, orient_key, item)
        c.execute(insert_string)


def get_name(header, splitchar="_", items=2, regex=None, repl=None):
    """parse the name of a locus from a file"""
    name = "_".join(header.split(splitchar)[:items]).lstrip('>').strip().lower()
    if regex is not None and repl is not None:
        return re.sub(regex, repl, name)
    else:
        return name


def get_dupes(lastz_file, regex=None, repl=None):
    """Given a lastz_file of probes aligned to themselves, get duplicates"""
    matches = defaultdict(list)
    dupes = set()
    for lz in lastz.Reader(lastz_file):
        target_name = get_name(lz.name1, "|", 1)
        query_name = get_name(lz.name2, "|", 1)
        matches[target_name].append(query_name)
    # see if one probe matches any other probes
    # other than the children of the locus
    for k, v in matches.iteritems():
        # if the probe doesn't match itself, we have
        # problems
        if len(v) > 1:
            for i in v:
                if i != k:
                    dupes.add(k)
                    dupes.add(i)
        elif k != v[0]:
            dupes.add(k)
    if not regex:
        return dupes
    else:
        return set([re.sub(regex, repl, d).lower() for d in dupes])


def contig_count(contig):
    """Return a count of contigs from a fasta file"""
    return sum([1 for line in open(contig, 'rU').readlines() if line.startswith('>')])


def get_organism_names_from_fasta_files(ff):
    """Given a fasta file name, parse taxon name from file name"""
    return [os.path.basename(f).split('.')[0].replace('-', "_")
        for f in ff]


def check_contigs_for_dupes(matches):
    """check for contigs that match more than 1 UCE locus"""
    node_dupes = defaultdict(list)
    for node in matches:
        node_dupes[node] = len(set(matches[node]))
    dupe_set = set([node for node in node_dupes if node_dupes[node] > 1])
    return dupe_set


def check_probes_for_dupes(revmatches):
    """Check for UCE probes that match more than one contig"""
    dupe_set = set([i for uce, node in revmatches.iteritems()
                if len(node) > 1 for i in list(node)])
    return dupe_set


def pretty_print_output(critter, matches, contigs, pd, mc, mp):
    """Write some nice output to stdout"""
    unique_matches = sum([1 for node, uce in matches.iteritems()])
    out = "\t {0}: {1} ({2:.2f}%) uniques of {3} contigs, {4} dupe probe matches, " + \
            "{5} UCE probes matching multiple contigs, {6} contigs matching multiple UCE probes"
    print out.format(
            critter,
            unique_matches,
            float(unique_matches) / contigs * 100,
            contigs,
            len(pd),
            len(mp),
            len(mc)
        )


def get_contig_name(header):
    """parse the contig name from the header of either velvet/trinity assembled contigs"""
    match = re.search("^(Node_\d+|comp\d+_c\d+_seq\d+).*", header)
    return match.groups()[0]


def main():
    args = get_args()
    if args.regex and args.repl is not None:
        # "s_[0-9]+$"
        regex = re.compile(args.regex)
        uces = set([get_name(read.identifier, "|", 1, regex=regex, repl=args.repl)
            for read in fasta.FastaReader(args.query)])
    else:
        uces = set([get_name(read.identifier, "|", 1)
            for read in fasta.FastaReader(args.query)])
        regex = None
    if args.dupefile:
        print "\t Getting dupes"
        dupes = get_dupes(args.dupefile, regex, args.repl)
    fasta_files = glob.glob(os.path.join(args.contigs, '*.fa*'))
    organisms = get_organism_names_from_fasta_files(fasta_files)
    conn, c = create_probe_database(
            os.path.join(args.output, 'probe.matches.sqlite'),
            organisms,
            uces
        )
    print "Processing:"
    for contig in fasta_files:
        critter = os.path.basename(contig).split('.')[0].replace('-', "_")
        output = os.path.join(
                    args.output, \
                    os.path.splitext(os.path.basename(contig))[0] + '.lastz'
                )
        contigs = contig_count(contig)
        # align the probes to the contigs
        alignment = lastz.Align(
                contig,
                args.query,
                args.coverage,
                args.identity,
                output
            )
        lzstdout, lztstderr = alignment.run()
        # parse the lastz results of the alignment
        matches, orientation, revmatches = \
                defaultdict(set), defaultdict(set), defaultdict(set)
        probe_dupes = set()
        if not lztstderr:
            for lz in lastz.Reader(output):
                # get strandedness of match
                contig_name = get_contig_name(lz.name1)
                uce_name = get_name(lz.name2, "|", 1, regex=regex, repl=args.repl)
                if args.dupefile and uce_name in dupes:
                    probe_dupes.add(uce_name)
                else:
                    matches[contig_name].add(uce_name)
                    orientation[uce_name].add(lz.strand2)
                    revmatches[uce_name].add(contig_name)
        # we need to check nodes for dupe matches to the same probes
        contigs_matching_mult_uces = check_contigs_for_dupes(matches)
        uces_matching_mult_contigs = check_probes_for_dupes(revmatches)
        nodes_to_drop = contigs_matching_mult_uces.union(uces_matching_mult_contigs)
        # remove dupe and/or dubious nodes/contigs
        match_copy = copy.deepcopy(matches)
        for k in match_copy.keys():
            if k in nodes_to_drop:
                del matches[k]
        store_lastz_results_in_db(c, matches, orientation, critter)
        conn.commit()
        pretty_print_output(
                critter,
                matches,
                contigs,
                probe_dupes,
                contigs_matching_mult_uces,
                uces_matching_mult_contigs
            )

if __name__ == '__main__':
    main()
