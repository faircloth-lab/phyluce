#!/usr/bin/env python
# encoding: utf-8

"""
get_align_summary_data.py

Created by Brant Faircloth on 16 August 2011.
Copyright 2011 Brant C. Faircloth. All rights reserved.

The program iterates through a folder of nexus files and returns the count
of alignments having more than "--percent" of "--taxa" taxa.
"""


import os
import math
import numpy
import argparse
import multiprocessing
from Bio import AlignIO
from collections import Counter

from phyluce.helpers import is_dir, FullPaths, get_alignment_files
from phyluce.log import setup_logging

import pdb


def get_args():
    parser = argparse.ArgumentParser(
        description="""Compute summary statistics for alignments in parallel""",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "--alignments",
        required=True,
        type=is_dir,
        action=FullPaths,
        help="""The directory containing alignments to be summarized."""
    )
    parser.add_argument(
        "--input-format",
        dest="input_format",
        choices=["fasta", "nexus", "phylip", "phylip-relaxed", "clustal", "emboss", "stockholm"],
        default="nexus",
        help="""The input alignment format.""",
    )
    parser.add_argument(
        "--show-taxon-counts",
        action="store_true",
        default=False,
        help="""Show the count of loci with X taxa.""",
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
        "--cores",
        type=int,
        default=1,
        help="""Process alignments in parallel using --cores for alignment. """ +
        """This is the number of PHYSICAL CPUs."""
    )
    return parser.parse_args()


class AlignMeta:
    def __init__(self):
        self.length = None
        self.taxa = None
        self.missing = None
        self.gaps = None
        self.characters = None
        self.nucleotides = None


def get_characters(aln, nucleotides):
    cnt = Counter()
    percent_missing = []
    for seq in aln:
        seq_string = str(seq.seq).upper()
        cnt.update(seq_string)
        percent_missing.append(float(seq_string.count("?")) / len(seq_string))
    return cnt, numpy.mean(numpy.array(percent_missing))


def get_stats(work):
    file, format = work
    aln = AlignIO.read(file, format)
    nucleotides = set(["A", "C", "G", "T"])
    meta = AlignMeta()
    meta.length = aln.get_alignment_length()
    meta.taxa = len(aln)
    meta.characters, meta.percent_missing = get_characters(aln, nucleotides)
    meta.nucleotides = Counter()
    for k, v in meta.characters.iteritems():
        if k in nucleotides:
            meta.nucleotides.update({k:v})
    meta.gaps = meta.characters["-"]
    meta.missing = meta.characters["?"]
    return meta


def get_lengths(summary):
    lengths = numpy.array([aln.length for aln in summary])
    total = numpy.sum(lengths)
    mean = numpy.mean(lengths)
    ci  = 1.96 * (numpy.std(lengths, ddof=1) / numpy.sqrt(len(lengths)))
    min = numpy.min(lengths)
    max = numpy.max(lengths)
    return total, mean, ci, min, max


def get_taxa(summary):
    taxa = numpy.array([aln.taxa for aln in summary])
    mean = numpy.mean(taxa)
    ci  = 1.96 * (numpy.std(taxa, ddof=1) / numpy.sqrt(len(taxa)))
    min = numpy.min(taxa)
    max = numpy.max(taxa)
    cnt = Counter(taxa)
    return cnt, mean, ci, min, max


def get_percent_missing(summary):
    missing = numpy.array([aln.percent_missing for aln in summary])
    mean = numpy.mean(missing)
    ci  = 1.96 * (numpy.std(missing, ddof=1) / numpy.sqrt(len(missing)))
    min = numpy.min(missing)
    max = numpy.max(missing)
    return mean, ci, min, max


def total_characters(summary):
    all = Counter()
    for aln in summary:
        all.update(aln.characters)
    return all, sum(all.values())


def total_nucleotides(summary):
    all = Counter()
    for aln in summary:
        all.update(aln.nucleotides)
    return sum(all.values())


def get_matrix_percentages(t_cnt):
    # get max taxa in alignments
    mx = max(t_cnt.keys())
    # get percentages
    stops = {}
    for i in numpy.arange(0.5, 1, 0.05):
        stops[i] = math.floor(i * mx)
    percentages = {}
    for percent, stop in stops.iteritems():
        total = 0
        for cnt, aln in t_cnt.iteritems():
            if cnt >= stop:
                total += aln
        percentages[percent] = total
    return percentages


def log_length_summary(log, loci, a_vars):
    a_total, a_mean, a_ci, a_min, a_max = a_vars
    text = " Alignment summary "
    log.info(text.center(65, "-"))
    log.info("[Alignments] loci:\t{:,}".format(loci))
    log.info("[Alignments] length:\t{:,}".format(a_total))
    log.info("[Alignments] mean:\t{:.2f}".format(a_mean))
    log.info("[Alignments] 95% CI:\t{:.2f}".format(a_ci))
    log.info("[Alignments] min:\t{}".format(a_min))
    log.info("[Alignments] max:\t{:,}".format(a_max))


def log_taxa_summary(log, t_vars):
    t_cnt, t_mean, t_ci, t_min, t_max = t_vars
    text = " Taxon summary "
    log.info(text.center(65, "-"))
    log.info("[Taxa] mean:\t\t{:.2f}".format(t_mean))
    log.info("[Taxa] 95% CI:\t{:.2f}".format(t_ci))
    log.info("[Taxa] min:\t\t{}".format(t_min))
    log.info("[Taxa] max:\t\t{}".format(t_max))


def log_missing_summary(log, m_vars):
    m_mean, m_ci, m_min, m_max = m_vars
    text = " Missing data from trim summary "
    log.info(text.center(65, "-"))
    log.info("[Missing] mean:\t{:.2f}".format(m_mean * 100))
    log.info("[Missing] 95% CI:\t{:.2f}".format(m_ci * 100))
    log.info("[Missing] min:\t{:.2f}".format(m_min * 100))
    log.info("[Missing] max:\t{:.2f}".format(m_max * 100))


def log_char_summary(log, sum_characters, sum_nucleotides):
    text = " Character count summary "
    log.info(text.center(65, "-"))
    log.info("[All characters]\t{:,}".format(sum_characters))
    log.info("[Nucleotides]\t\t{:,}".format(sum_nucleotides))


def log_matrix_summary(log, percentages):
    text = " Data matrix completeness summary "
    log.info(text.center(65, "-"))
    for k in sorted(percentages.keys()):
        log.info("[Matrix {0}%]\t\t{1} alignments".format(
            int(k * 100),
            percentages[k],
        ))


def log_taxa_dist(log, show_taxon_counts, t_cnt):
    if show_taxon_counts:
        text = " Alignment counts by taxa present "
        log.info(text.center(65, "-"))
        for k in sorted(t_cnt.keys()):
            log.info("[Taxa] {0} alignments contain {1:,} taxa".format(
                t_cnt[k],
                k,
            ))


def log_character_dist(log, all_bases):
    text = " Character counts "
    log.info(text.center(65, "-"))
    for k in sorted(all_bases.keys()):
        if k in ['A','C','G','T','a','c','g','t', '-', '?']:
            log.info("[Characters] '{0}' is present {1:,} times".format(
                k,
                all_bases[k],
            ))
        else:
            log.warn("[Characters] '{0}' is present {1:,} times".format(
                k,
                all_bases[k],
            ))


def main():
    args = get_args()
    # setup logging
    log, my_name = setup_logging(args)
    # find all alignments
    files = get_alignment_files(log, args.alignments, args.input_format)
    work = [[file, args.input_format] for file in files]
    log.info("Computing summary statistics using {} cores".format(args.cores))
    if args.cores > 1:
        assert args.cores <= multiprocessing.cpu_count(), "You've specified more cores than you have"
        pool = multiprocessing.Pool(args.cores)
        summary = pool.map(get_stats, work)
    else:
        summary = map(get_stats, work)
    # alignments
    a_vars = get_lengths(summary)
    log_length_summary(log, len(summary), a_vars)
    # taxa
    t_vars = get_taxa(summary)
    log_taxa_summary(log, t_vars)
    # missing
    m_vars = get_percent_missing(summary)
    log_missing_summary(log, m_vars)
    # characters
    all_bases, sum_characters = total_characters(summary)
    sum_nucleotides = total_nucleotides(summary)
    log_char_summary(log, sum_characters, sum_nucleotides)
    # matrix
    percentages = get_matrix_percentages(t_vars[0])
    log_matrix_summary(log, percentages)
    # taxa dist.
    log_taxa_dist(log, args.show_taxon_counts, t_vars[0])
    # character dist
    log_character_dist(log, all_bases)
    # end
    text = " Completed {} ".format(my_name)
    log.info(text.center(65, "="))


if __name__ == '__main__':
    main()
