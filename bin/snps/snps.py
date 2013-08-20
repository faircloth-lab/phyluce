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


def create_bwa_index_files(log, reference):
    log.info("Running bwa indexing against {}".format(reference))
    cwd = os.getcwd()
    # move into reference directory
    os.chdir(os.path.dirname(reference))
    cmd = ["bwa", "index", reference]
    with open('bwa-index-file.txt', 'a') as outf:
        proc = subprocess.Popen(cmd, stdout=outf, stderr=subprocess.STDOUT)
        proc.communicate()
    # mvoe back to working directory
    os.chdir(cwd)


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
            create_bwa_index_files(log, reference)
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


def run_create_sai(log, sample, sample_dir, ref, cores, reads, read):
    log.info("Creating read index file for {}".format(reads.file))
    cmd = [
        "bwa",
        "aln",
        "-t",
        str(cores),
        ref,
        reads.pth
    ]
    aln_out_fname = os.path.join(sample_dir, '{}-r{}.sai'.format(sample, read))
    aln_err_fname = os.path.join(sample_dir, '{}-r{}.bwa-aln-out.txt'.format(sample, read))
    with open(aln_out_fname, 'w') as aln_out:
        with open(aln_err_fname, 'w') as aln_err:
            proc = subprocess.Popen(cmd, stdout=aln_out, stderr=aln_err)
            proc.communicate()
    return aln_out_fname


def bwa_create_sai(log, sample, sample_dir, ref, cores, reads, read):
    if read == 1:
        aln_out_fname = run_create_sai(log, sample, sample_dir, ref, cores, reads, read)
    elif read == 2:
        aln_out_fname = run_create_sai(log, sample, sample_dir, ref, cores, reads, read)
    return aln_out_fname


def align_read1_and_read2_with_bwa(log, sample, sample_dir, ref, cores, r1, r2):
    r1sai = bwa_create_sai(log, sample, sample_dir, ref, cores, r1, 1)
    r2sai = bwa_create_sai(log, sample, sample_dir, ref, cores, r2, 2)
    cmd1 = [
        "bwa",
        "sampe",
        "-a",
        "700",
        ref,
        r1sai,
        r2sai,
        r1.pth,
        r2.pth
    ]
    cmd2 = [
        "samtools",
        "view",
        "-bS",
        "-"
    ]
    sampe_out_fname = os.path.join(sample_dir, '{}.bwa-sampe-out.txt'.format(sample))
    samtools_out_fname = os.path.join(sample_dir, '{}.samtools-out.txt'.format(sample))
    bam_out_fname = os.path.join(sample_dir, '{}.bam'.format(sample))
    log.info("Building BAM for {}".format(sample))
    with open(sampe_out_fname, 'w') as sampe_out:
        with open(samtools_out_fname, 'w') as samtools_out:
            with open(bam_out_fname, 'w') as bam_out:
                proc1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE, stderr=sampe_out)
                proc2 = subprocess.Popen(cmd2, stdin=proc1.stdout, stdout=bam_out, stderr=samtools_out)
                proc1.stdout.close()
                proc2.communicate()
    # remove the sai files (they'll be stale soon)
    os.remove(r1sai)
    os.remove(r2sai)
    return bam_out_fname


def new_bam_name(bam, append):
    pth, bamfname = os.path.split(bam)
    bamfname = os.path.splitext(bamfname)[0]
    new_bamfname = "{}-{}.bam".format(bamfname, append)
    new_bam = os.path.join(pth, new_bamfname)
    return new_bam


def clean_up_bam(log, sample, sample_dir, bam):
    log.info("Cleaning BAM for {}".format(sample))
    new_bam = new_bam_name(bam, "CL")
    cmd = [
        "java",
        "-Xmx20g",
        "-jar",
        "/home/bcf/bin/CleanSam.jar",
        "I={}".format(bam),
        "O={}".format(new_bam)
    ]
    picard_clean_out_fname = os.path.join(sample_dir, '{}.picard-clean-out.txt'.format(sample))
    with open(picard_clean_out_fname, 'w') as picard_out:
        proc = subprocess.Popen(cmd, stdout=picard_out, stderr=subprocess.STDOUT)
        proc.communicate()
    # remove old bam
    os.remove(bam)
    return new_bam

def add_rg_header_info(log, sample, sample_dir, flowcell, bam):
    #pdb.set_trace()
    log.info("Adding RG header to BAM for {}".format(sample))
    new_bam = new_bam_name(bam, "RG")
    cmd = [
        "java",
        "-Xmx20g",
        "-jar",
        "/home/bcf/bin/AddOrReplaceReadGroups.jar",
        "I={}".format(bam),
        "O={}".format(new_bam),
        "SORT_ORDER=coordinate",
        "RGPL=illumina",
        "RGPU={}".format(flowcell),
        "RGLB=Lib1",
        "RGID={}".format(sample),
        "RGSM={}".format(sample),
        "VALIDATION_STRINGENCY=LENIENT"
    ]
    picard_rg_out_fname = os.path.join(sample_dir, '{}.picard-RG-out.txt'.format(sample))
    with open(picard_rg_out_fname, 'w') as picard_out:
        proc = subprocess.Popen(cmd, stdout=picard_out, stderr=subprocess.STDOUT)
        proc.communicate()
    # remove old bam
    os.remove(bam)
    return new_bam


def mark_and_remove_dupes(log, sample, sample_dir, bam):
    log.info("Removing read duplicates from BAM for {}".format(sample))
    new_bam = new_bam_name(bam, "DD")
    metricsfile = os.path.join(sample_dir, "{}.picard-metricsfile.txt".format(sample))
    cmd = [
        "java",
        "-Xmx20g",
        "-jar",
        "/home/bcf/bin/MarkDuplicates.jar",
        "I={}".format(bam),
        "O={}".format(new_bam),
        "METRICS_FILE={}".format(metricsfile),
        "MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=250",
        "ASSUME_SORTED=true",
        "VALIDATION_STRINGENCY=SILENT",
        "REMOVE_DUPLICATES=true",
    ]
    picard_dd_out_fname = os.path.join(sample_dir, '{}.picard-DD-out.txt'.format(sample))
    with open(picard_dd_out_fname, 'w') as picard_out:
        proc = subprocess.Popen(cmd, stdout=picard_out, stderr=subprocess.STDOUT)
        proc.communicate()
    # remove old bam
    os.remove(bam)
    return new_bam


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
            bam = align_read1_and_read2_with_bwa(log, sample, sample_dir, reference, args.cores, fastq.r1, fastq.r2)
            # clean the bam up (MAPq 0 and trim overlapping reads)
            bam = clean_up_bam(log, sample, sample_dir, bam)
            # get flowcell id
            fc = flowcells[sample]
            bam = add_rg_header_info(log, sample, sample_dir, fc, bam)
            if not args.no_remove_duplicates:
                bam = mark_and_remove_dupes(log, sample, sample_dir, bam)
            else:
                log.info("You have selected to keep apparent duplicate reads")
        #if fastq.singleton:
            # bwa align singleton reads
        #    pass

if __name__ == '__main__':
    main()
