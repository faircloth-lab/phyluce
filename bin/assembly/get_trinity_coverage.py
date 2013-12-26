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
import glob
import shutil
import argparse
from phyluce.helpers import FullPaths, is_file, is_dir
from phyluce.third_party import which
from phyluce.raw_reads import get_input_data, get_fastq_input_files
from phyluce.bwa import *
from phyluce.log import setup_logging

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
        "--log-path",
        action=FullPaths,
        type=is_dir,
        default=None,
        help="""The path to a directory to hold logs."""
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
    parser.add_argument(
        "--velvet",
        action="store_true",
        default=False,
        help="""Compute coverage across velvet/abYss contigs instead of Trinity contigs""",
    )
    return parser.parse_args()


def cleanup_trinity_assembly_folder(log, pth):
    log.info("Removing extraneous Trinity files")
    files = glob.glob(os.path.join(pth, '*'))
    for file in files:
        if not os.path.basename(file) in ("Trinity.fasta", "trinity.log"):
            if os.path.isfile(file):
                os.remove(file)
            elif os.path.isdir(file):
                shutil.rmtree(file)


def symlink_trimmed_contigs(log, sample, contig_dir, trimmed_fasta_path):
    log.info("Symlinking trimmed contigs into contigs-trimmed")
    try:
        linkpth = os.path.join(contig_dir, "{}.contigs.fasta".format(sample))
        relpth = os.path.relpath(trimmed_fasta_path, linkpth)
        os.symlink(relpth, linkpth)
    except:
        log.warn("Unable to symlink {} to {}".format(sample, linkpth))


def main():
    # get args and options
    args = get_args()
    # setup logging
    log, my_name = setup_logging(args)
    # get the input data
    log.info("Getting input filenames")
    input = get_input_data(args.assemblo_config, None)
    # Get path to bwa
    try:
        bwa = which('bwa')[0]
    except:
        raise EnvironmentError("Cannot find bwa.  Ensure it is installed and in your $PATH")
    # make the symlink directory within the output directory
    contig_dir = os.path.join(args.assemblies, 'contigs-trimmed')
    if not os.path.isdir(contig_dir):
        os.makedirs(contig_dir)
    else:
        pass
    for group in input:
        sample, reads = group
        # pretty print taxon status
        text = " Processing {} ".format(sample)
        log.info(text.center(65, "-"))
        # ensure that assembly exists
        assembly_pth = os.path.join(args.assemblies, sample)
        assembly = os.path.join(assembly_pth, "contigs.fasta")
        if not os.path.exists(assembly):
            raise IOError("Assembly for {} does not appear to exist.".format(sample))
        if args.clean:
            cleanup_trinity_assembly_folder(log, assembly_pth)
        # determine the types of raw read data that we have
        fastq = get_fastq_input_files(reads, args.subfolder, log)
        # create the bwa index
        bwa_create_index_files(log, assembly)
        samtools_create_faidx(log, sample, assembly_pth, assembly)
        picard_create_reference_dict(log, sample, assembly_pth, assembly)
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
        # get singleton reads for alignment
        if args.bwa_mem and fastq.singleton:
            bam_se = bwa_mem_se_align(log, sample, assembly_pth, assembly, args.cores, fastq.singleton)
            bam_se = picard_clean_up_bam(log, sample, assembly_pth, bam_se, 'se')
            bam_se = picard_add_rg_header_info(log, sample, assembly_pth, "Generic", bam_se, "se")
        # if we only have se reads, those will be in fastq.r1 only
        elif args.bwa_mem and not fastq.r2 and fastq.r1:
            bam_se = bwa_mem_se_align(log, sample, assembly_pth, assembly, args.cores, fastq.r1)
            bam_se = picard_clean_up_bam(log, sample, assembly_pth, bam_se, 'se')
            bam_se = picard_add_rg_header_info(log, sample, assembly_pth, "Generic", bam_se, "se")
        elif not args.bwa_mem and fastq.singleton:
            bam_se = bwa_se_align(log, sample, assembly_pth, assembly, args.cores, fastq.singleton)
            bam_se = picard_clean_up_bam(log, sample, assembly_pth, bam_se, 'se')
            bam_se = picard_add_rg_header_info(log, sample, assembly_pth, "Generic", bam_se, "se")
        elif not args.bwa_mem and not fastq.r2 and fastq.r1:
            bam_se = bwa_se_align(log, sample, assembly_pth, assembly, args.cores, fastq.r1)
            bam_se = picard_clean_up_bam(log, sample, assembly_pth, bam_se, 'se')
            bam_se = picard_add_rg_header_info(log, sample, assembly_pth, "Generic", bam_se, "se")
        if bam and bam_se:
            bam = picard_merge_two_bams(log, sample, assembly_pth, bam, bam_se)
        elif bam_se and not bam:
            bam = bam_se
        if not bam:
            raise IOError("There is no BAM file.  Check bwa log files for problems.")
        samtools_index(log, sample, assembly_pth, bam)
        coverage = gatk_coverage(log, sample, assembly_pth, assembly, args.cores, bam)
        overall_contigs = get_coverage_from_gatk(log, sample, assembly_pth, coverage, args.velvet)
        remove_gatk_coverage_files(log, assembly_pth, coverage)
        trimmed_fasta_path = filter_screened_contigs_from_assembly(log, sample, assembly_pth, assembly, overall_contigs)
        symlink_trimmed_contigs(log, sample, contig_dir, trimmed_fasta_path)
    # end
    text = " Completed {} ".format(my_name)
    log.info(text.center(65, "="))




if __name__ == '__main__':
    main()
