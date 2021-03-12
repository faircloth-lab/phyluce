#!/usr/bin/env python
# encoding: utf-8

"""
helpers.py

Created by Brant Faircloth on 15 July 2011.
Copyright 2011 Brant C. Faircloth. All rights reserved.
"""

import os
import re
import sys
import glob
import argparse
import shutil
import configparser
from collections import defaultdict

from rich import prompt

from phyluce import lastz
from phyluce.pth import get_all_user_params

# import pdb


class FullPaths(argparse.Action):
    """Expand user- and relative-paths"""

    def __call__(self, parser, namespace, values, option_string=None):
        setattr(
            namespace, self.dest, os.path.abspath(os.path.expanduser(values))
        )


class CreateDir(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        # get the full path
        d = os.path.abspath(os.path.expanduser(values))
        # check to see if directory exists
        if os.path.exists(d):
            if prompt.Confirm.ask(
                "[magenta][WARNING] Output directory exists, REMOVE[/magenta]"
            ):
                shutil.rmtree(d)
            else:
                print("[QUIT]")
                sys.exit()
        # create the new directory
        os.makedirs(d)
        # return the full path
        setattr(namespace, self.dest, d)


class CreateFile(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        # get the full path
        f = os.path.abspath(os.path.expanduser(values))
        # check to see if directory exists
        if os.path.exists(f):
            if prompt.Confirm.ask(
                "[magenta][WARNING] Output file exists, REMOVE[/magenta]"
            ):
                os.remove(f)
            else:
                print("[QUIT]")
                sys.exit()
        setattr(namespace, self.dest, f)


def is_dir(dirname):
    if not os.path.isdir(dirname):
        msg = "{0} is not a directory".format(dirname)
        raise argparse.ArgumentTypeError(msg)
    else:
        return dirname


def is_file(filename):
    if not os.path.isfile:
        msg = "{0} is not a file".format(filename)
        raise argparse.ArgumentTypeError(msg)
    else:
        return filename


def get_name(header, splitchar="_", items=2):
    """use own function vs. import from match_contigs_to_probes - we don't want lowercase"""
    if splitchar:
        return "_".join(header.split(splitchar)[:items]).lstrip(">")
    else:
        return header.lstrip(">")


def get_dupe_matches(lastz_file, splitchar="|", pos=1, longfile=False):
    matches = defaultdict(list)
    for lz in lastz.Reader(lastz_file, longfile):
        target_name = get_name(lz.name1, splitchar, pos)
        query_name = get_name(lz.name2, splitchar, pos)
        matches[target_name].append(query_name)
    return matches


def get_dupes(lastz_file, splitchar="|", pos=1, longfile=False):
    dupes = set()
    matches = get_dupe_matches(lastz_file, splitchar, pos, longfile)
    # see if one probe matches any other probes
    # other than the children of the locus
    for k, v in matches.items():
        # if the probe doesn't match itself, we have
        # problems
        if len(v) > 1:
            for i in v:
                if i != k:
                    dupes.add(k)
                    dupes.add(i)
        elif k != v[0]:
            dupes.add(k)
    return dupes


def get_names_from_config(config, group):
    try:
        return [i[0] for i in config.items(group)]
    except configparser.NoSectionError:
        return None


def get_file_extensions(ftype):
    ext = {
        "fasta": (".fasta", ".fsa", ".aln", ".fa"),
        "nexus": (".nexus", ".nex"),
        "phylip": (".phylip", ".phy"),
        "phylip-relaxed": (".phylip", ".phy", ".phylip-relaxed"),
        "phylip-sequential": (".phylip", ".phy", ".phylip-sequential"),
        "clustal": (".clustal", ".clw"),
        "emboss": (".emboss",),
        "stockholm": (".stockholm",),
    }
    return ext[ftype]


def get_alignment_files(log, input_dir, input_format):
    log.info("Getting alignment files")
    alignments = []
    for ftype in get_file_extensions(input_format):
        alignments.extend(
            glob.glob(os.path.join(input_dir, "*{}".format(ftype)))
        )
    if not alignments:
        log.critical(
            "No alignment files found.  Check --input-format argument."
        )
    return alignments


def write_alignments_to_outdir(log, outdir, alignments, output_format):
    log.info("Writing output files")
    for tup in alignments:
        locus, aln = tup
        if aln.trimmed is not None:
            outname = "{}{}".format(
                os.path.join(outdir, locus),
                get_file_extensions(output_format)[0],
            )
            with open(outname, "w") as outf:
                outf.write(format(aln.trimmed, output_format))
        else:
            log.warn("DROPPED {0} from output".format(locus))


def get_contig_header_string():
    return "|".join(get_all_user_params("headers"))
