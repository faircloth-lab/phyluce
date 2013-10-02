#!/usr/bin/env python
# encoding: utf-8
"""
File: get_trinity_coverage.py
Author: Brant Faircloth

Created by Brant Faircloth on 30 September 2013 10:09 PDT (-0700)
Copyright (c) 2013 Brant C. Faircloth. All rights reserved.

Description:

"""

import os
import sys
import glob
import shutil
import logging
import argparse
from phyluce.helpers import FullPaths, is_file
from phyluce.third_party import which
from phyluce.raw_reads import get_input_data, get_fastq_input_files
from phyluce.bwa import *

import pdb

def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
        description="""Get coverage of UCE assemblies from Trinity"""
    )
    parser.add_argument(
        "--assemblies",
        required=True,
        action=FullPaths,
        default=None,
        help="""The directory containing the assemblies"""
    )
    parser.add_argument(
        "--assemblo-config",
        type=is_file,
        action=FullPaths,
        default=None,
        help="""The assemblo_trinity configuration file"""
    )
    parser.add_argument(
        "--cores",
        type=int,
        default=1,
        help="""The number of compute cores/threads to run with Trinity"""
    )
    parser.add_argument(
        "--subfolder",
        type=str,
        default='',
        help="""A subdirectory, below the level of the group, containing the reads"""
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
        "--bwa-mem",
        action="store_true",
        default=False,
        help="""Use bwa-mem instead of standard bwa""",
    )
    return parser.parse_args()


def setup_logging(level):
    log = logging.getLogger("Trinity Coverage")
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


def cleanup_trinity_assembly_folder(log, pth):
    log.info("Removing extraneous Trinity files")
    files = glob.glob(os.path.join(pth, '*'))
    for file in files:
        if not os.path.basename(file) in ("Trinity.fasta", "trinity.log"):
            if os.path.isfile(file):
                os.remove(file)
            elif os.path.isdir(file):
                shutil.rmtree(file)

def main():
    # get args and options
    args = get_args()
    # setup logger
    log = setup_logging(args.verbosity)
    log.info("=================== Starting Trinity Coverage ===================")
    # get the input data
    log.info("Getting input filenames")
    input = get_input_data(args.assemblo_config, None)
    # Get path to bwa
    try:
        bwa = which('bwa')[0]
    except:
        raise EnvironmentError("Cannot find bwa.  Ensure it is installed and in your $PATH")
    for group in input:
        sample, reads = group
        # pretty print taxon status
        text = " Processing {} ".format(sample)
        log.info(text.center(65, "-"))
        # ensure that assembly exists
        assembly_pth = os.path.join(args.assemblies, sample)
        assembly = os.path.join(assembly_pth, "Trinity.fasta")
        if not os.path.exists(assembly):
            raise IOError("Assembly for {} does not appear to exist.".format(sample))
        if args.clean:
            cleanup_trinity_assembly_folder(log, assembly_pth)
        # determine the types of raw read data that we have
        fastq = get_fastq_input_files(reads, args.subfolder, log)
        # create the bwa index
        bwa_create_index_files(log, assembly)
        bam = False
        bam_se = False
        if args.bwa_mem and fastq.r1 and fastq.r2:
            bam = bwa_mem_pe_align(log, sample, assembly_pth, assembly, args.cores, fastq.r1, fastq.r2)
            bam = picard_clean_up_bam(log, sample, assembly_pth, bam, "pe")
            bam = picard_add_rg_header_info(log, sample, assembly_pth, "Generic", bam, "pe")
        elif not args.bwa_mem and fastq.r1 and fastq.r2:
            bam = bwa_pe_align(log, sample, assembly_pth, assembly, args.cores, fastq.r1, fastq.r2)
            bam = picard_clean_up_bam(log, sample, assembly_pth, bam, "pe")
            bam = picard_add_rg_header_info(log, sample, assembly_pth, "Generic", bam, "pe")
        if args.bwa_mem and fastq.singleton:
            bam_se = bwa_mem_se_align(log, sample, assembly_pth, assembly, args.cores, fastq.singleton)
            bam_se = picard_clean_up_bam(log, sample, assembly_pth, bam_se, 'se')
            bam_se = picard_add_rg_header_info(log, sample, assembly_pth, "Generic", bam_se, "se")
        elif not args.bwa_mem and fastq.singleton:
            bam_se = bwa_se_align(log, sample, assembly_pth, assembly, args.cores, fastq.singleton)
            bam_se = picard_clean_up_bam(log, sample, assembly_pth, bam_se, 'se')
            bam_se = picard_add_rg_header_info(log, sample, assembly_pth, "Generic", bam_se, "se")
        if bam and bam_se:
            bam = picard_merge_two_bams(log, sample, assembly_pth, bam, bam_se)
        elif bam_se and not bam:
            bam = bam_se
        samtools_index(log, sample, assembly_pth, bam)
        coverage = gatk_coverage(log, sample, assembly_pth, assembly, args.cores, bam)
        overall_contigs = get_coverage_from_gatk(log, sample, assembly_pth, coverage)
        remove_gatk_coverage_files(log, assembly_pth, coverage)
        filter_screened_contigs_from_assembly(log, sample, assembly_pth, assembly, overall_contigs)




if __name__ == '__main__':
    main()
