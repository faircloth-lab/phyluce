#!/usr/bin/env python
# encoding: utf-8
"""
File: slice_sequence_from_genomes2.py
Author: Brant Faircloth

Created by Brant Faircloth on 08 June 2013 11:06 PDT (-0700)
Copyright (c) 2013 Brant C. Faircloth. All rights reserved.

Description:

"""

import os
import re
import sys
import copy
import logging
import argparse
import ConfigParser
from collections import defaultdict

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from bx.seq import twobit
from phyluce import lastz
from phyluce.helpers import FullPaths, is_dir, is_file


import pdb

def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
        description="""Given a LASTZ input directory, find matches, add flank, and return a FASTA file of sequences""")
    parser.add_argument(
        "conf",
        action=FullPaths,
        type=is_file,
        help="""Path to the configuration file"""
    )
    parser.add_argument(
        "lastz",
        action=FullPaths,
        type=is_dir,
        help="""Path to the directory containing LASTZ results"""
    )
    parser.add_argument(
        "output",
        action=FullPaths,
        help="""Path to the output directory for storing FASTA files"""
    )
    parser.add_argument(
        "--flank",
        type=int,
        default=500,
        help="""The amount of flanking sequence to add to each match""",
    )
    parser.add_argument(
        "--name-pattern",
        dest="pattern",
        type=str,
        default=None,
        help="An alternate name pattern to transform the conf entry into"
    )
    parser.add_argument(
        "--exclude",
        type=str,
        nargs='+',
        default=None,
        help="""Species to exclude from genome slicing""",
    )
    parser.add_argument(
        "--verbosity",
        type=str,
        choices=["INFO", "WARN", "CRITICAL"],
        default="INFO",
        help="""The logging level to use"""
    )
    parser.add_argument(
        "--contig_orient",
        action="store_true",
        default=False,
        help="""Check orientation by contigs versus probes - useful for multi-species probe sets""",
    )
    return parser.parse_args()


def get_all_files_from_conf(conf, pattern=None):
    all_files = []
    if conf.has_section("chromos"):
        all_files.extend(conf.items("chromos"))
    if conf.has_section("scaffolds"):
        all_files.extend(conf.items("scaffolds"))
    if pattern is not None:
        files = [(v[0], pattern.format(v[0]), v[1]) for v in all_files]
    else:
        files = [(v[0], v[0], v[1]) for v in all_files]
    return files


def new_get_probe_name(header, regex="^(uce-\d+)(?:_p\d+.*)"):
    match = re.search(regex, header)
    return match.groups()[0]


def check_loci_for_dupes(matches):
    """Check for UCE loci that match more than one contig"""
    dupe_set = set([uce for uce, contigs in matches.iteritems() if len(contigs) > 1])
    return dupe_set


def slice_and_return_fasta(tb, name, min, max, flank):
    if min - flank > 0:
        ss = min - flank
    else:
        ss = 0
    if max + flank < len(tb[name]):
        se = max + flank
    else:
        se = len(tb[name])
    return ss, se, tb[name][ss:se]


def remove_ambiguous_ends(ss, se, sequence):
    start_matches = re.search("^([Nn]+)", sequence)
    end_matches = re.search("([Nn]+)$", sequence)
    if start_matches:
        new_start = len(start_matches.groups()[0])
        ss = ss + new_start
    else:
        new_start = 0
    if end_matches:
        new_end = len(end_matches.groups()[0])
        se = se - new_end
    else:
        new_end = 0
    return ss, se, sequence[new_start:len(sequence) - new_end]


def remove_repetitive_ends(ss, se, sequence):
    start_matches = re.search("^([acgt]+)", sequence)
    end_matches = re.search("([acgt]+)$", sequence)
    if start_matches:
        new_start = len(start_matches.groups()[0])
        ss = ss + new_start
    else:
        new_start = 0
    if end_matches:
        new_end = len(end_matches.groups()[0])
        se = se - new_end
    else:
        new_end = 0
    return ss, se, sequence[new_start:len(sequence) - new_end]


def build_sequence_object(cnt, contig, ss, se, uce, min, max, orient, sorted_positions, sequence):
    ss, se, sequence = remove_ambiguous_ends(ss, se, sequence)
    ss, se, sequence = remove_repetitive_ends(ss, se, sequence)
    name = "Node_{0}_length_{1}_cov_1000 |contig:{2}|slice:{3}-{4}|uce:{5}|match:{6}-{7}|orient:{8}|probes:{9}".format(
        cnt,
        len(sequence),
        contig,
        ss,
        se,
        uce,
        min,
        max,
        list(orient)[0],
        len(sorted_positions)
    )
    return SeqRecord(Seq(sequence), id=name, name='', description='')


def parse_lastz_file(lz, contig_orient):
    all_uce_names = set()
    uce_matches = defaultdict(lambda: defaultdict(list))
    orientation = defaultdict(lambda: defaultdict(set))
    for lz in lastz.Reader(lz, long_format=True):
        # get strandedness of match
        contig_name = lz.name1
        # get name of UCE from lastz info
        uce_name = new_get_probe_name(lz.name2)
        # keep a record of all UCEs matched
        all_uce_names.add(uce_name)
        uce_matches[uce_name][contig_name].append([lz.zstart1, lz.end1])
        # usually, we get orientation matches by probes. but for mixed species
        # probes, they may be in several orientations, which can cause problems,
        # so check their orientation relative to the contig
        if contig_orient:
            orientation[uce_name][contig_name].add(lz.strand1)
        else:
            orientation[uce_name][contig_name].add(lz.strand2)
    return all_uce_names, uce_matches, orientation


def setup_logging(level):
    log = logging.getLogger("Phyluce")
    console = logging.StreamHandler(sys.stdout)
    if level == "INFO":
        log.setLevel(logging.INFO)
        console.setLevel(logging.INFO)
    if level == "WARN":
        log.setLevel(logging.WARN)
        console.setLevel(logging.WARN)
    if level == "CRITICAL":
        log.setLevel(logging.CRITICAL)
        console.setLevel(logging.CRITICAL)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    console.setFormatter(formatter)
    log.addHandler(console)
    return log


def main():
    args = get_args()
    # make the output directory
    try:
        os.makedirs(args.output)
    except OSError:
        answer = raw_input("Output directory exists.  Overwrite [Y/n]? ")
        if answer != 'Y':
            sys.exit()
    # setup logger
    log = setup_logging(args.verbosity)
    log.info("=================== Starting Phyluce: Slice Sequence ===================")
    # parse config file of genome locations (genomes in 2bit format)
    conf = ConfigParser.ConfigParser()
    conf.optionxform = str
    conf.read(args.conf)
    # get the files associated with the config entries
    all_files = get_all_files_from_conf(conf, args.pattern)
    for genome in all_files:
        short_name, long_name, twobit_name = genome
        text = " Working on {} genome ".format(short_name)
        log.info(text.center(65, "-"))
        if not args.exclude or (short_name not in args.exclude):
            out_name = os.path.join(args.output, "{}.fasta".format(short_name.lower()))
            with open(out_name, 'w') as outf:
                log.info("Reading {} genome".format(short_name))
                tb = twobit.TwoBitFile(file(twobit_name))
                # parse the lastz results of the alignment
                lz = os.path.join(args.lastz, long_name)
                all_uce_names, uce_matches, orientation = parse_lastz_file(lz, args.contig_orient)
                # we need to check nodes for dupe matches to the same probes
                uce_loci_matching_mult_contigs = check_loci_for_dupes(uce_matches)
                # delete those loci that hit multiple probes or had multiple probes hit them
                # this does not filter out contigs hitting multiple loci for the simple
                # reason that many loci will hit the same "contig" in circumstances where
                # the contig is large/chromosome sized.  These will get filtered out in
                # the next step of the matching process after we've created fasta sequences
                # representing each UCE locus
                for k in uce_matches.keys():
                    if k in uce_loci_matching_mult_contigs:
                        del uce_matches[k]
                # QC contig matches
                node_count = 0
                orient_drop = set()
                length_drop = set()
                for uce_name, matches in uce_matches.iteritems():
                    bad = False
                    # make sure there is only one UCE match
                    assert len(matches.keys()) == 1, "There are multiple UCE matches"
                    for contig_name, positions in matches.iteritems():
                        # remove any probes with mixed orientation
                        orient = orientation[uce_name][contig_name]
                        if len(orient) > 1:
                            bad = True
                            orient_drop.add(uce_name)
                        if not bad:
                            # sort the positions for each contig
                            sorted_positions = sorted(positions)
                            if len(sorted_positions) > 1:
                                for i in range(1, len(sorted_positions)):
                                    # drop those contigs where probes fall > 500 bp apart
                                    if sorted_positions[i][0] - sorted_positions[i-1][1] > 500:
                                        bad = True
                                        length_drop.add(uce_name)
                                        break
                                    else:
                                        bad = False
                        if not bad:
                            min = sorted_positions[0][0]
                            max = sorted_positions[-1][-1]
                            ss, se, sequence = slice_and_return_fasta(tb, contig_name, min, max, args.flank)
                            seq = build_sequence_object(node_count, contig_name, ss, se, uce_name, min, max, orient, sorted_positions, sequence)
                            outf.write(seq.format('fasta'))
                            node_count += 1
            output = "{}: {} uces, {} dupes, {} non-dupes, {} orient drop, {} length drop, {} written".format(
                short_name,
                len(all_uce_names),
                len(uce_loci_matching_mult_contigs),
                len(uce_matches.keys()),
                len(orient_drop),
                len(length_drop),
                node_count
            )
            log.info(output)


if __name__ == '__main__':
    main()
