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
from phyluce.third_party import which
from phyluce.helpers import FullPaths, is_dir, is_file
from phyluce.assembly import get_fastq_input_files
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
        help="""The number of compute cores/threads to run with Trinity"""
    )
    parser.add_argument(
        "--verbosity",
        type=str,
        choices=["INFO", "WARN", "CRITICAL"],
        default="INFO",
        help="""The logging level to use"""
    )
    parser.add_argument(
        "--clean",
        action="store_true",
        default=False,
        help="""Cleanup all intermediate Trinity files""",
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


def setup_logging(level):
    log = logging.getLogger("Assemblo")
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
    # get args and options
    args = get_args()
    # setup logger
    log = setup_logging(args.verbosity)
    text = " Starting BWA alignment "
    log.info(text.center(85, "-"))
    # get the config file data
    conf = ConfigParser.ConfigParser(allow_no_value=True)
    conf.optionxform = str
    conf.read(args.config)
    # get the input data
    log.info("Getting input filenames and creating output directories")
    reference, individuals = get_input_data(log, conf, args.output)
    flowcells = dict(conf.items("flowcell"))
    for indiv in individuals:
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
            bam = picard_clean_up_bam(log, sample, sample_dir, bam)
            # get flowcell id
            fc = flowcells[sample]
            bam = picard_add_rg_header_info(log, sample, sample_dir, fc, bam)
            if not args.no_remove_duplicates:
                bam = picard_mark_and_remove_dupes(log, sample, sample_dir, bam)
            else:
                log.info("You have selected to keep apparent duplicate reads")
        #if fastq.singleton:
            # bwa align singleton reads
        #    pass

if __name__ == '__main__':
    main()
