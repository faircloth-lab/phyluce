# -*- coding: utf-8 -*-

"""
(c) 2015 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 22 May 2015 14:56 CDT (-0500)
"""

import os
import math
import numpy
from Bio import AlignIO
from collections import Counter
from phyluce import sites


class AlignMeta:
    def __init__(self):
        self.length = None
        self.taxa = None
        self.missing = None
        self.gaps = None
        self.characters = None
        self.nucleotides = None
        self.sum_informative_sites = None
        self.sum_differences = None
        self.sum_counted_sites = None


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
    meta.name = os.path.basename(file)
    meta.length = aln.get_alignment_length()
    meta.taxa = len(aln)
    meta.characters, meta.percent_missing = get_characters(aln, nucleotides)
    meta.nucleotides = Counter()
    for k, v in meta.characters.iteritems():
        if k in nucleotides:
            meta.nucleotides.update({k:v})
    meta.gaps = meta.characters["-"]
    meta.missing = meta.characters["?"]
    meta.sum_informative_sites, meta.sum_differences, meta.sum_counted_sites = sites.compute_informative_sites(aln)
    return meta


def get_lengths(summary):
    lengths = numpy.array([aln.length for aln in summary])
    total = numpy.sum(lengths)
    mean = numpy.mean(lengths)
    ci  = 1.96 * (numpy.std(lengths, ddof=1) / numpy.sqrt(len(lengths)))
    min = numpy.min(lengths)
    max = numpy.max(lengths)
    return total, mean, ci, min, max


def get_sites(summary):
    sites = numpy.array([aln.sum_informative_sites for aln in summary])
    total = numpy.sum(sites)
    mean = numpy.mean(sites)
    ci  = 1.96 * (numpy.std(sites, ddof=1) / numpy.sqrt(len(sites)))
    min = numpy.min(sites)
    max = numpy.max(sites)
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
        # add a little fudge factor to deal with floats
        stops[i] = math.ceil((i - 0.01) * mx)
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


def log_sites_summary(log, loci, s_vars):
    s_total, s_mean, s_ci, s_min, s_max = s_vars
    text = " Informative Sites summary "
    log.info(text.center(65, "-"))
    log.info("[Sites] loci:\t{:,}".format(loci))
    log.info("[Sites] total:\t{:,}".format(s_total))
    log.info("[Sites] mean:\t{:.2f}".format(s_mean))
    log.info("[Sites] 95% CI:\t{:.2f}".format(s_ci))
    log.info("[Sites] min:\t{}".format(s_min))
    log.info("[Sites] max:\t{:,}".format(s_max))


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
