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
from phyluce.helpers import is_dir, is_file, FullPaths
from phyluce.log import setup_logging
from collections import defaultdict
from Bio import SeqIO

#import pdb


def get_args():
    parser = argparse.ArgumentParser(description="Match UCE probes/baits to assembled contigs " +
        "and store the data in a relational database.  The matching process is dependent on the " +
        "probe names in the file.  If the probe names are not like 'uce-1001_p1' where 'uce-' " +
        "indicates we're searching for uce loci, '1001' indicates locus 1001, '_p1' indicates " +
        "this is probe 1 for locus 1001, you will need to set the optional --regex parameter. " +
        "So, if your probe names are 'MyProbe-A_probe1', the --regex will look like " +
        "--regex='^(MyProbe-\W+)(?:_probe\d+.*)'"
    )
    parser.add_argument(
        '--contigs',
        required=True,
        type=is_dir,
        action=FullPaths,
        help="The directory containing the assembled contigs in which you are searching for UCE loci."
    )
    parser.add_argument(
        '--probes',
        required=True,
        type=is_file,
        action=FullPaths,
        help="The bait/probe file in FASTA format."
    )
    parser.add_argument(
        '--output',
        required=True,
        action=FullPaths,
        help="The directory in which to store the resulting SQL database and LASTZ files."
    )
    parser.add_argument(
        "--verbosity",
        type=str,
        choices=["INFO", "WARN", "CRITICAL"],
        default="INFO",
        help="""The logging level to use."""
    )
    parser.add_argument(
        "--log-path",
        action=FullPaths,
        type=is_dir,
        default=None,
        help="""The path to a directory to hold logs."""
    )
    parser.add_argument(
        '--min-coverage',
        default=80,
        type=int,
        help="The minimum percent coverage required for a match [default=80]."
    )
    parser.add_argument(
        '--min-identity',
        default=80,
        type=int,
        help="The minimum percent identity required for a match [default=80]."
    )
    parser.add_argument(
        '--dupefile',
        help="Path to self-to-self lastz results for baits to remove potential duplicate probes."
    )
    parser.add_argument(
        "--regex",
        type=str,
        default="^(uce-\d+)(?:_p\d+.*)",
        help="""A regular expression to apply to the probe names for replacement [default='^(uce-\d+)(?:_p\d+.*)'].""",
    )
    parser.add_argument(
        "--keep-duplicates",
        type=str,
        default=None,
        help="""A optional output file in which to store those loci that appear to be duplicates.""",
    )
    args = parser.parse_args()
    return args


def create_probe_database(log, db, organisms, uces):
    """Create the UCE-match database"""
    log.info("Creating the UCE-match database")
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
        log.critical("Database already exists")
        if e[0] == 'table matches already exists':
            answer = raw_input("Database already exists.  Overwrite [Y/n]? ")
            if answer == "Y" or "YES":
                os.remove(db)
                conn, c = create_probe_database(db, organisms, uces)
            else:
                sys.exit(2)
        else:
            log.critical("Cannot create database")
            raise sqlite3.OperationalError("Cannot create database")
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


def get_dupes(log, lastz_file, regex):
    """Given a lastz_file of probes aligned to themselves, get duplicates"""
    log.info("Checking probe/bait sequences for duplicates")
    matches = defaultdict(list)
    dupes = set()
    # get names and strip probe designation since loci are the same
    for lz in lastz.Reader(lastz_file):
        target_name = new_get_probe_name(lz.name1, regex)
        query_name = new_get_probe_name(lz.name2, regex)
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
    # make sure all names are lowercase
    return set([d.lower() for d in dupes])


def contig_count(contig):
    """Return a count of contigs from a fasta file"""
    return sum([1 for line in open(contig, 'rU').readlines() if line.startswith('>')])


def get_organism_names_from_fasta_files(ff):
    """Given a fasta file name, parse taxon name from file name"""
    return [os.path.basename(f).split('.')[0].replace('-', "_") for f in ff]


def check_contigs_for_dupes(matches):
    """check for contigs that match more than 1 UCE locus"""
    node_dupes = defaultdict(list)
    for node in matches:
        node_dupes[node] = len(set(matches[node]))
    dupe_set = set([node for node in node_dupes if node_dupes[node] > 1])
    return dupe_set


def check_loci_for_dupes(revmatches):
    """Check for UCE probes that match more than one contig"""
    dupe_contigs = []
    dupe_uces = []
    for uce, node in revmatches.iteritems():
        if len(node) > 1:
            dupe_contigs.extend(node)
            dupe_uces.append(uce)
    #dupe_contigs = set([i for uce, node in revmatches.iteritems() if len(node) > 1 for i in list(node)])
    #pdb.set_trace()
    return set(dupe_contigs), set(dupe_uces)


def pretty_log_output(log, critter, matches, contigs, pd, mc, uce_dupe_uces):
    """Write some nice output to the logfile/stdout"""
    unique_matches = sum([1 for node, uce in matches.iteritems()])
    out = "{0}: {1} ({2:.2f}%) uniques of {3} contigs, {4} dupe probe matches, " + \
        "{5} UCE loci removed for matching multiple contigs, {6} contigs " + \
        "removed for matching multiple UCE loci"
    log.info(
        out.format(
            critter,
            unique_matches,
            float(unique_matches) / contigs * 100,
            contigs,
            len(pd),
            len(uce_dupe_uces),
            len(mc)
        )
    )


def get_contig_name(header):
    """parse the contig name from the header of either velvet/trinity assembled contigs"""
    match = re.search("^(node_\d+|comp\d+_c\d+_seq\d+|c\d+_g\d+_i\d+).*", header, flags=re.I)
    return match.groups()[0]


def new_get_probe_name(header, regex):
    match = re.search(regex, header)
    return match.groups()[0]


def main():
    args = get_args()
    log, my_name = setup_logging(args)
    regex = re.compile(args.regex)
    if not os.path.isdir(args.output):
        os.makedirs(args.output)
    else:
        raise IOError("The directory {} already exists.  Please check and remove by hand.".format(args.output))
    uces = set(new_get_probe_name(seq.id, regex) for seq in SeqIO.parse(open(args.probes, 'rU'), 'fasta'))
    if args.dupefile:
        dupes = get_dupes(log, args.dupefile, regex)
    else:
        dupes = set()
    fasta_files = glob.glob(os.path.join(args.contigs, '*.fa*'))
    organisms = get_organism_names_from_fasta_files(fasta_files)
    conn, c = create_probe_database(
        log,
        os.path.join(args.output, 'probe.matches.sqlite'),
        organisms,
        uces
    )
    log.info("Processing contig data")
    # open a file for duplicate writing, if we're interested
    if args.keep_duplicates is not None:
        dupefile = open(args.keep_duplicates, 'w')
    else:
        dupefile = None
    log.info("{}".format("-" * 65))
    for contig in sorted(fasta_files):
        critter = os.path.basename(contig).split('.')[0].replace('-', "_")
        output = os.path.join(
            args.output,
            os.path.splitext(os.path.basename(contig))[0] + '.lastz'
        )
        contigs = contig_count(contig)
        # align the probes to the contigs
        alignment = lastz.Align(
            contig,
            args.probes,
            args.min_coverage,
            args.min_identity,
            output
        )
        lzstdout, lztstderr = alignment.run()
        if lztstderr:
            raise EnvironmentError("lastz: {}".format(lztstderr))
        # parse the lastz results of the alignment
        matches = defaultdict(set)
        orientation = defaultdict(set)
        revmatches = defaultdict(set)
        probe_dupes = set()
        if not lztstderr:
            for lz in lastz.Reader(output):
                # get strandedness of match
                contig_name = get_contig_name(lz.name1)
                uce_name = new_get_probe_name(lz.name2, regex)
                if args.dupefile and uce_name in dupes:
                    probe_dupes.add(uce_name)
                else:
                    matches[contig_name].add(uce_name)
                    orientation[uce_name].add(lz.strand2)
                    revmatches[uce_name].add(contig_name)
        # we need to check nodes for dupe matches to the same probes
        contigs_matching_mult_uces = check_contigs_for_dupes(matches)
        uce_dupe_contigs, uce_dupe_uces = check_loci_for_dupes(revmatches)
        nodes_to_drop = contigs_matching_mult_uces.union(uce_dupe_contigs)
        # write out duplicates if requested
        if dupefile is not None:
            log.info("Writing duplicates file for {}".format(critter))
            if len(uce_dupe_uces) != 0:
                dupefile.write("[{} - probes hitting multiple contigs]\n".format(critter))
                for uce in uce_dupe_uces:
                    dupefile.write("{}:{}\n".format(uce, ', '.join(revmatches[uce])))
                dupefile.write("\n")
            if len(contigs_matching_mult_uces) != 0:
                dupefile.write("[{} - contigs hitting multiple probes]\n".format(critter))
                for dupe in contigs_matching_mult_uces:
                    dupefile.write("{}:{}\n".format(dupe, ', '.join(matches[dupe])))
                dupefile.write("\n")
        #pdb.set_trace()
        # remove dupe and/or dubious nodes/contigs
        match_copy = copy.deepcopy(matches)
        for k in match_copy.keys():
            if k in nodes_to_drop:
                del matches[k]
        store_lastz_results_in_db(c, matches, orientation, critter)
        conn.commit()
        pretty_log_output(
            log,
            critter,
            matches,
            contigs,
            probe_dupes,
            contigs_matching_mult_uces,
            uce_dupe_uces
        )
    if dupefile is not None:
        dupefile.close()
    log.info("{}".format("-" * 65))
    log.info("The LASTZ alignments are in {}".format(args.output))
    log.info("The UCE match database is in {}".format(os.path.join(args.output, "probes.matches.sqlite")))
    text = " Completed {} ".format(my_name)
    log.info(text.center(65, "="))

if __name__ == '__main__':
    main()
