#!/usr/bin/env python
# encoding: utf-8
"""
File: assemble_contigs_from_bgi_genomes.py
Author: Brant Faircloth

Created by Brant Faircloth on 29 March 2012 15:03 PDT (-0700)
Copyright (c) 2012 Brant C. Faircloth. All rights reserved.

Description: 

"""

import os
import re
import sys
import glob
import sqlite3
import tempfile
import argparse

from phyluce import lastz
from operator import itemgetter
from phyluce.muscle import Align
from collections import defaultdict
from seqtools.sequence import fasta
from phyluce.helpers import FullPaths, is_dir, is_file, get_dupes, get_name, run_checks, snip_if_many_N_bases, get_uce_names_from_probes


import pdb


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
            description="""Assemble fake (velvet) contigs from matches to faircloth+stephen probes from BGI""")
    parser.add_argument('probes',
            action=FullPaths,
            type=is_file,
            help='Fasta file of the probe sequence searched',
        )
    parser.add_argument('lastz',
            action=FullPaths,
            type=is_dir,
            help='Location of the LASTZ results',
        )
    parser.add_argument(
            "fasta",
            action=FullPaths,
            type=is_dir,
            help="""Location of the input fasta slices"""
        )
    parser.add_argument(
            "output",
            action=FullPaths,
            help="""Location for the output fake contigs"""
        )
    parser.add_argument(
            "--xref",
            action=FullPaths,
            help="""Location of xref mapping file from BGI"""
        )
    parser.add_argument('--flank',
            help='The flank length to add',
            default=500,
            type=int
        )
    parser.add_argument(
            "--db",
            action=FullPaths,
            type=is_file,
            default=None,
            help="""Add columns to existing database rather than starting over""",
        )
    parser.add_argument(
            "--extend",
            action="store_true",
            default=False,
            help="""Add columns to existing database rather than starting  over""",
        )
    parser.add_argument(
            "--oldprobe",
            action="store_true",
            default=False,
            help="""Generating a database for old probe naming scheme""",
        )
    parser.add_argument('--dupefile',
            action=FullPaths,
            type=is_file,
            help='The path to a lastz file of lastz-against-self results'
        )
    return parser.parse_args()


def get_taxon_from_filename(filename):
    return os.path.basename(filename).split('.')[0].lower()


def get_taxa_names_from_fastas(fastas):
    taxa = []
    for ff in glob.glob(os.path.join(fastas, '*')):
        taxon = get_taxon_from_filename(ff)
        taxa.append(taxon)
    return taxa


def create_bgi_name_map(xref, regex):
    probes = {}
    locus = {}
    for line in open(xref, 'rU'):
        ls = line.strip().split('\t')
        probes[ls[1]] = ls[0]
        locus[re.sub(regex, '', ls[1])] = '|'.join(ls[0].split('|')[:2])
    return probes, locus


def get_fasta_name_from_lastz_pth(lastz, fasta):
    name = os.path.splitext(os.path.basename(lastz))[0]
    try:
        for ext in ['.fa', '.fasta', '.gz', '.fasta.gz', '.fa.gz']:
            tempname = os.path.join(fasta, os.path.splitext(name)[0]) + ext
            if os.path.isfile(tempname):
                break
    except:
        raise ValueError("Cannot find the fasta file {}".format(name))
    return tempname


def get_bgi_matches(lastz_file, stripnum):
    matches = defaultdict(list)
    probes = defaultdict(int)
    for lz in lastz.Reader(lastz_file, long_format=True):
        uce_name = re.sub(stripnum, 's', lz.name2).lower()
        probe_number = int(lz.name2.split('_')[-1])
        if probe_number > probes[uce_name]:
            probes[uce_name] = probe_number
        matches[uce_name].append(
                [
                    get_name(lz.name1).lower(),
                    lz.strand2,
                    lz.zstart1,
                    lz.end1
                ]
            )
    return matches, probes


def get_old_probe_matches(lastz_file):
    matches = defaultdict(list)
    probes = defaultdict(int)
    for lz in lastz.Reader(lastz_file, long_format=True):
        uce_name = lz.name2.split('|')[0]
        probe_number = int(lz.name2.split(':')[-1])
        if probe_number > probes[uce_name]:
            probes[uce_name] = probe_number
        matches[uce_name].append(
                [
                    get_name(lz.name1).lower(),
                    lz.strand2,
                    lz.zstart1,
                    lz.end1
                ]
            )
    return matches, probes


def quality_control_matches(matches, probes, dupes, k, v, verbose=False):
    """check to make sure we don't get any more matches than expected
    and that matches are reasonably close to each other on their respective
    chromos"""
    skip = []
    if k in dupes:
        skip.append(k)
        if verbose:
            print "{0} is in dupefile".format(k)
    elif len(v) > 1:
        if run_checks(k, v, probes, verbose=False):
            # sort by match position
            v_sort = sorted(v, key=itemgetter(2))
            start, end = v_sort[0][2], v_sort[-1][3]
            diff = end - start
            # ensure our range is less than N(probes) * probe_length - this
            # still gives us a little wiggle room because probes are ~ 2X tiled
            if diff > (probes[k] * 140):
                skip.append(k)
                if verbose:
                    print "range longer than expected"
        else:
            skip.append(k)
    return skip


def get_probe_positions(record):
    """given a fasta record, get probe position from info"""
    idsplit = record.identifier.split('|')
    ps, pe = [int(i) for i in idsplit[2].split('-')]
    rs, re = [int(i) for i in idsplit[3].split('-')]
    probe_start = ps - rs
    probe_end = pe - rs
    return probe_start, probe_end


def trim_uce_reads(record, flank):
    """given a fasta record, slice from record based on probe position
    and flanking sequence length"""
    idsplit = record.identifier.split('|')
    ps, pe = [int(i) for i in idsplit[2].split('-')]
    rs, re = [int(i) for i in idsplit[3].split('-')]
    probe_start = ps - rs
    probe_end = pe - rs
    if probe_start - flank > 0:
        new_start = probe_start - flank
    else:
        new_start = 0
    if probe_end + flank < re - rs:
        new_end = probe_end + flank
    else:
        new_end = re - rs
    new_read_start = rs + new_start
    new_read_end = rs + new_end
    record.identifier = "{0}|{1}|{2}-{3}|{4}-{5}||{6}|{7}|{8}".format(
            idsplit[0],
            idsplit[1],
            ps,
            pe,
            new_read_start,
            new_read_end,
            idsplit[5],
            idsplit[6],
            idsplit[7]
        )
    record.sequence = record.sequence[new_start:new_end]
    return record


def create_probe_database(db, organisms, uces, genome=False):
    """docstring for create_probe_database"""
    conn = sqlite3.connect(db)
    c = conn.cursor()
    c.execute("PRAGMA foreign_keys = ON")
    try:
        create_string = [org + ' text' for org in organisms]
        query = "CREATE TABLE matches (uce text primary key, {0})".format(','.join(create_string))
        c.execute(query)
        query = "CREATE TABLE match_map (uce text primary key, {0})".format(','.join(create_string))
        c.execute(query)
        if genome:
            query = "CREATE TABLE contig_map (species text, locus text, contig text, new_name text)"
            c.execute(query)
        for uce in uces:
            c.execute("INSERT INTO matches(uce) values (?)", (uce,))
            c.execute("INSERT INTO match_map(uce) values (?)", (uce,))
    except sqlite3.OperationalError, e:
        if e[0] == 'table matches already exists':
            answer = raw_input("Database already exists.  Overwrite [Y/n]? ")
            if answer == "Y" or "YES":
                os.remove(db)
                conn, c = create_probe_database(db, organisms, uces, genome)
            else:
                sys.exit(2)
        else:
            raise sqlite3.OperationalError
            pdb.set_trace()
    return conn, c


def extend_probe_database(db, organisms):
    conn = sqlite3.connect(db)
    c = conn.cursor()
    c.execute("PRAGMA foreign_keys = ON")
    org_string = [org + ' text' for org in organisms]
    for org in org_string:
        for table in ["matches", "match_map"]:
            try:
                query = "ALTER TABLE {0} ADD COLUMN {1}".format(table, org)
                c.execute(query)
            except sqlite3.OperationalError, e:
                if "duplicate column name" in e[0]:
                    print "Resetting {} for {}".format(table, org)
                    query = "UPDATE {0} SET {1} = NULL".format(table, org.strip(" text"))
                    c.execute(query)
    return conn, c


def main():
    args = get_args()
    # compile some regular expressions we'll use later
    stripnum = re.compile("s_[0-9]+$")
    manyn = re.compile("[N,n]{20,}")
    # get names of loci and taxa
    uces = get_uce_names_from_probes(args.probes, regex=stripnum, repl='s', lower=True)
    taxa = get_taxa_names_from_fastas(args.fasta)
    print "\n"
    if not args.extend:
        if args.db is None:
            db = os.path.join(args.output, 'probe.matches.sqlite')
        else:
            db = args.db
        # create db to hold results
        conn, c = create_probe_database(
                db,
                taxa,
                uces,
                True
            )
    else:
        conn, c = extend_probe_database(
                args.db,
                taxa
            )
    # get duplicate probe sequences for filtering
    if args.dupefile:
        print "Determining duplicate probes..."
        dupes = get_dupes(args.dupefile, longfile=False)
    else:
        dupes = None
    # because of structure, strip probe designation from dupes
    # leaving only locus name.  lowercase all.
    dupes = set([re.sub(stripnum, 's', d).lower() for d in dupes])
    # iterate over LASTZ files for each taxon
    for lz in glob.glob(os.path.join(args.lastz, '*')):
        # get taxon name from lastz file
        taxon = get_taxon_from_filename(lz)
        print "\n{0}\n{1}\n{0}".format('=' * 30, taxon)
        # get fasta name from lastz file
        ff = get_fasta_name_from_lastz_pth(lz, args.fasta)
        # get lastz matches
        print "\tGetting LASTZ matches from GENOME alignments..."
        if not args.oldprobe:
            matches, probes = get_bgi_matches(lz, stripnum)
        else:
            matches, probes = get_old_probe_matches(lz)
        # remove bad loci (dupes)
        print "\tGetting bad (potentially duplicate) GENOME matches..."
        loci_to_skip = []
        for k, v in matches.iteritems():
            # check matches to makes sure all is well - keep names lc
            loci_to_skip.extend(quality_control_matches(matches, probes, dupes, k, v, False))
        #pdb.set_trace()
        # convert to set, to keep only uniques
        loci_to_skip = set(loci_to_skip)
        print "\tSkipping {} bad (duplicate hit) loci...".format(len(loci_to_skip))
        # get (and possibly assemble) non-skipped
        seqdict = defaultdict(list)
        # determine those contigs to skip and group those to assemble
        for contig in fasta.FastaReader(ff):
            # make sure all names are lowercase
            contig.identifier = contig.identifier.lower()
            if not args.oldprobe:
                name = contig.identifier.split('|')[-3]
                locus = re.sub(stripnum, 's', name)
            else:
                locus = contig.identifier.split('|')[-5]
            # skip what we identified as bad loci
            if locus not in loci_to_skip:
                seqdict[locus].append(contig)
        output_name = "{}.fasta".format(taxon.replace('_', '-'))
        fout_name = os.path.join(args.output, output_name)
        print "\tOutput filename is {}".format(output_name)
        fout = fasta.FastaWriter(fout_name)
        # this tracks "fake" contig number
        count = 0
        # this tracks loci kept
        kept = 0
        # when > 1 contig, assemble contigs across matches
        sys.stdout.write("\tWriting and Aligning/Assembling UCE loci with multiple probes (dot/1000 loci)")
        for k, v in seqdict.iteritems():
            bad = False
            contig_names = []
            if count % 1000 == 0:
                sys.stdout.write('.')
                sys.stdout.flush()
            if len(v) == 1:
                # trim ambiguous bases on flanks
                record = v[0]
                orient = [matches[k][0][1]]
                if args.flank:
                    record = trim_uce_reads(record, args.flank)
                contig_names.append(record.identifier)
                record.sequence = record.sequence.strip('N')
                # trim many ambiguous bases within contig
                result = manyn.search(record.sequence)
                if result:
                    uce_start, uce_end = get_probe_positions(record)
                    uce = record.sequence[uce_start:uce_end]
                    record.sequence = snip_if_many_N_bases(manyn, k, record.sequence, uce, verbose=False)
                # change header
                record.identifier = ">Node_{0}_length_{1}_cov_1000".format(
                        count,
                        len(record.sequence)
                    )
                fout.write(v[0])
            else:
                orient = list(set([m[1] for m in matches[k]]))
                # skip any loci having matches of mixed orientation
                # ['+', '-']
                if len(orient) == 1:
                    # create tempfile for the reads
                    fd, temp = tempfile.mkstemp(suffix='.fasta')
                    os.close(fd)
                    temp_out = fasta.FastaWriter(temp)
                    # write all slices to outfile, trimming if we want
                    #pdb.set_trace()
                    for record in v:
                        if args.flank:
                            record = trim_uce_reads(record, args.flank)
                        # keep names of contigs we assembled to store in db assoc
                        # w/ resulting assembled contig name
                        contig_names.append(record.identifier)
                        record.sequence = record.sequence.strip('N')
                        # trim many ambiguous bases within contig
                        result = manyn.search(record.sequence)
                        if result:
                            uce_start, uce_end = get_probe_positions(record)
                            uce = record.sequence[uce_start:uce_end]
                            record.sequence = snip_if_many_N_bases(manyn, k, record.sequence, uce, verbose=False)
                        temp_out.write(record)
                    # make sure to close the file
                    temp_out.close()
                    # assemble
                    aln = Align(temp)
                    aln.run_alignment()
                    record = fasta.FastaSequence()
                    record.sequence = aln.alignment_consensus.tostring()
                    record.identifier = ">Node_{0}_length_{1}_cov_1000".format(
                            count,
                            len(record.sequence)
                        )
                    fout.write(record)
                else:
                    bad = True
            if not bad:
                # track contig assembly and renaming data in db
                q = "UPDATE matches SET {0} = 1 WHERE uce = '{1}'".format(taxon, k)
                c.execute(q)
                # generate db match and match map tables for data
                orient_key = "node_{0}({1})".format(count, orient[0])
                q = "UPDATE match_map SET {0} = '{1}' WHERE uce = '{2}'".format(taxon, orient_key, k)
                c.execute(q)
                # keep track of new name :: old name mapping
                for old_name in contig_names:
                    q = "INSERT INTO contig_map VALUES ('{0}', '{1}', '{2}', '{3}')".format(taxon, k, old_name, record.identifier)
                    c.execute(q)
                kept += 1
            # tracking "fake" contig number
            count += 1
        conn.commit()
        print "\n\t{0} loci of {1} matched ({2:.0f}%), {3} dupes dropped ({4:.0f}%), {5} ({6:.0f}%) kept".format(
            count,
            len(uces),
            float(count) / len(uces) * 100,
            len(loci_to_skip),
            float(len(loci_to_skip)) / len(uces) * 100,
            kept,
            float(kept) / len(uces) * 100
            )
    #conn.commit()
    c.close()
    conn.close()

if __name__ == '__main__':
    main()
