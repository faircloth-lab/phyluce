#!/usr/bin/env python
# encoding: utf-8
"""
File: snps.py
Author: Brant Faircloth

Created by Brant Faircloth on 29 July 2013 16:07 PDT (-0700)
Copyright (c) 2013 Brant C. Faircloth. All rights reserved.

Description:

"""


import os
import sys
import logging
import argparse
import subprocess
import ConfigParser
from phyluce.log import setup_logging
from phyluce.third_party import which
from phyluce.helpers import FullPaths, is_dir, is_file
from phyluce.raw_reads import get_fastq_input_files
from phyluce.bwa import *

import pdb


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
        description="""Call SNPs"""
    )
    parser.add_argument(
        "--config",
        required=True,
        type=is_file,
        action=FullPaths,
        default=None,
        help="""A configuration file containing"""
    )
    parser.add_argument(
        "--output",
        required=True,
        action=FullPaths,
        default=None,
        help="""The directory in which to store the SNPs files"""
    )
    parser.add_argument(
        "--subfolder",
        type=str,
        default='',
        help="""A subdirectory, below the level of the group, containing the reads"""
    )
    parser.add_argument(
        "--cores",
        type=int,
        default=1,
        help="""The number of compute cores/threads to use"""
    )
    parser.add_argument(
        "--verbosity",
        type=str,
        choices=["INFO", "WARN", "CRITICAL"],
        default="INFO",
        help="""The logging level to use"""
    )
    parser.add_argument(
        "--log-path",
        action=FullPaths,
        type=is_dir,
        default=None,
        help="""The path to a directory to hold logs."""
    )
    parser.add_argument(
        "--no-remove-duplicates",
        action="store_true",
        default=False,
        help="""Do not remove duplicate reads.""",
    )
    return parser.parse_args()


def get_input_data(log, conf, output):
    # get reference sequence
    reference = conf.items('reference')
    # ensure there is 1 reference and it is a file
    assert len(reference) == 1, "There is more than one reference sequence listed."
    reference = reference[0][0]
    try:
        assert os.path.isfile(reference)
    except:
        raise IOError("{} is not a file".format(reference))
    # check reference to ensure that bwa has indexed
    for suffix in ['amb', 'ann', 'bwt', 'pac',  'sa']:
        bwa_file = "{}.{}".format(reference, suffix)
        try:
            assert os.path.isfile(bwa_file)
        except:
            log.info("Need to create BWA index file for reference")
            bwa_create_index_files(log, reference)
    individuals = conf.items('individuals')
    for sample in individuals:
        try:
            assert os.path.isdir(sample[1])
        except:
            raise IOError("{} is not a directory".format(sample[1]))
    return reference, individuals


def main():
    # get args and options
    args = get_args()
    # setup logging
    log, my_name = setup_logging(args)
    text = " Starting {} ".format(my_name)
    log.info(text.center(65, "="))
    # get the config file data
    conf = ConfigParser.ConfigParser(allow_no_value=True)
    conf.optionxform = str
    conf.read(args.config)
    # get the input data
    log.info("Getting input filenames and creating output directories")
    reference, individuals = get_input_data(log, conf, args.output)
    flowcells = dict(conf.items("flowcell"))
    for indiv in individuals:
        bam, bam_se = False, False
        sample, dir = indiv
        # pretty print taxon status
        text = " Processing {} ".format(sample)
        log.info(text.center(65, "-"))
        # make a directory for sample-specific assemblies
        sample_dir = os.path.join(args.output, sample)
        os.makedirs(sample_dir)
        # determine how many files we're dealing with
        fastq = get_fastq_input_files(dir, args.subfolder, log)
        if fastq.r1 and fastq.r2:
            # bwa align r1 and r2
            bam = bwa_pe_align(log, sample, sample_dir, reference, args.cores, fastq.r1, fastq.r2)
            # clean the bam up (MAPq 0 and trim overlapping reads)
            bam = picard_clean_up_bam(log, sample, sample_dir, bam, "pe")
            # get flowcell id
            fc = flowcells[sample]
            bam = picard_add_rg_header_info(log, sample, sample_dir, fc, bam, "pe")
            if not args.no_remove_duplicates:
                bam = picard_mark_and_remove_dupes(log, sample, sample_dir, bam, "pe")
            else:
                log.info("You have selected to keep apparent duplicate reads")
        if fastq.singleton:
            # bwa align singleton reads
            bam_se = bwa_se_align(log, sample, sample_dir, reference, args.cores, fastq.singleton)
            # clean the bam up (MAPq 0 and trim overlapping reads)
            bam_se = picard_clean_up_bam(log, sample, sample_dir, bam_se, "se")
            # get flowcell id
            fc = flowcells[sample]
            bam_se = picard_add_rg_header_info(log, sample, sample_dir, fc, bam_se, "se")
            if not args.no_remove_duplicates:
                bam_se = picard_mark_and_remove_dupes(log, sample, sample_dir, bam_se, "se")
            else:
                log.info("You have selected to keep apparent duplicate reads")
        if bam and bam_se:
            bam = picard_merge_two_bams(log, sample, sample_dir, bam, bam_se)
        elif bam_se and not bam:
            bam = bam_se
        if not bam:
            raise IOError("There is no BAM file.  Check bwa log files for problems.")
    # end
    text = " Completed {} ".format(my_name)
    log.info(text.center(65, "="))

if __name__ == '__main__':
    main()
