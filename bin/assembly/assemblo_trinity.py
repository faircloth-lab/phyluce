#!/usr/bin/env python
# encoding: utf-8
"""
File: assemblo_trinity.py
Author: Brant Faircloth

Created by Brant Faircloth on 02 June 2013 10:06 PDT (-0700)
Copyright (c) 2013 Brant C. Faircloth. All rights reserved.

Description: Assemble UCE cleaned read data using Trinity
(http://trinityrnaseq.sourceforge.net/).

Trinity is an assembly suite tailored to assembling data
where reads are at variable abundance.  This is in contrast
to other assemblers (e.g. Velvet, ABySS) that generally
assume a somewhat uniform coverage of reads across the genome
or "targets".

"""

import os
import glob
import shutil
import argparse
import subprocess
from phyluce.third_party import which
from phyluce.helpers import FullPaths, is_dir, is_file
from phyluce.raw_reads import get_input_data, get_input_files
from phyluce.log import setup_logging

#import pdb


def get_args():
    """Get arguments from CLI"""
    parser = argparse.ArgumentParser(
        description="""Assemble UCE raw read using Trinity"""
    )
    parser.add_argument(
        "--output",
        required=True,
        action=FullPaths,
        default=None,
        help="""The directory in which to store the assembly data"""
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
    # one of these is required.  The other will be set to None.
    input = parser.add_mutually_exclusive_group(required=True)
    input.add_argument(
        "--config",
        type=is_file,
        action=FullPaths,
        default=None,
        help="""A configuration file containing reads to assemble"""
    )
    input.add_argument(
        "--dir",
        type=is_dir,
        action=FullPaths,
        default=None,
        help="""A directory of reads to assemble""",
    )
    return parser.parse_args()


def copy_read_data(reads, sample_dir, log):
    """We need to combine singleton reads with read1
    data to include those data so trinity can use them"""
    # copy read data to sample dir
    log.info("Copying raw read data to {}".format(sample_dir))
    for fq in reads.reads:
        old_pth = os.path.join(fq.dir, fq.file)
        new_pth = os.path.join(sample_dir, fq.file)
        shutil.copy(old_pth, new_pth)
        # reset fq dirname
        fq.dir = sample_dir


def combine_read_data(reads, log):
    log.info("Combining singleton reads with R1 data")
    # setup cat command for subprocess
    cmd = [
        "cat",
        os.path.join(reads.singleton.dir, reads.singleton.file)
    ]
    # cat the singleton fastq/fasta contents into the r1 fastq/fasta file
    with open(os.path.join(reads.r1.dir, reads.r1.file), 'a') as outf:
        proc = subprocess.Popen(cmd, stdout=outf)
        proc.communicate()
    # remove the singleton file
    os.remove(os.path.join(reads.singleton.dir, reads.singleton.file))
    # unset the singleton fastq/fasta record
    reads.set_read('singleton', None, None)


def run_trinity_pe(trinity, reads, cores, log):
    log.info("Running Trinity.pl for PE data")
    cmd = [
        trinity,
        "--seqType",
        "fq",
        "--JM",
        "10G",
        "--left",
        os.path.join(reads.r1.dir, reads.r1.file),
        "--right",
        os.path.join(reads.r2.dir, reads.r2.file),
        "--CPU",
        str(cores),
        "--output",
        reads.r1.dir
    ]
    try:
        with open(os.path.join(reads.r1.dir, 'trinity.log'), 'w') as outf:
            proc = subprocess.Popen(cmd, stdout=outf)
            proc.communicate()
    except:
        log.critical("Could not assemble {}".format(reads.r1.dir))
    try:
        # trinity converts files to unzipped fastq/fasta. delete
        # gzips and fastas.  These data are also in `both.fa`
        # which is created at the start of Trinity.
        if reads.gzip:
            r1 = os.path.join(reads.r1.dir, reads.r1.file)
            r2 = os.path.join(reads.r2.dir, reads.r2.file)
            r1_fastq = os.path.splitext(r1)[0]
            r2_fastq = os.path.splitext(r2)[0]
            for f in [r1, r2, r1_fastq, r2_fastq]:
                os.remove(f)
        else:
            r1 = os.path.join(reads.r1.dir, reads.r1.file)
            r2 = os.path.join(reads.r2.dir, reads.r2.file)
            for f in [r1, r2]:
                os.remove(f)
    except:
        log.warn("Did not clean all fastq/fasta files from {}".format(reads.r1.dir))
    return reads.r1.dir


def run_trinity_se(trinity, reads, cores, log):
    log.info("Running Trinity.pl for SE data")
    cmd = [
        trinity,
        "--seqType",
        "fq",
        "--JM",
        "10G",
        "--single",
        os.path.join(reads.r1.dir, reads.r1.file),
        "--CPU",
        str(cores),
        "--output",
        reads.r1.dir
    ]
    try:
        with open(os.path.join(reads.r1.dir, 'trinity.log'), 'w') as outf:
            proc = subprocess.Popen(cmd, stdout=outf)
            proc.communicate()
    except:
        log.critical("Could not assemble {}".format(reads.r1.dir))
    try:
        # trinity converts files to unzipped fastq/fasta. delete
        # gzips and fastas.  These data are also in `both.fa`
        # which is created at the start of Trinity.
        if reads.gzip:
            r1 = os.path.join(reads.r1.dir, reads.r1.file)
            r1_fastq = os.path.splitext(r1)[0]
            for f in [r1, r1_fastq]:
                os.remove(f)
        else:
            r1 = os.path.join(reads.r1.dir, reads.r1.file)
            for f in [r1]:
                os.remove(f)
    except:
        log.warn("Did not clean all fastq/fasta files from {}".format(reads.r1.dir))
    return reads.r1.dir


def generate_symlinks(contig_dir, sample, reads, log):
    log.info("Symlinking assembled contigs into {}".format(contig_dir))
    try:
        trinity_fname = os.path.join(reads.r1.dir, "Trinity.fasta")
        contig_lname = os.path.join(reads.r1.dir, "contigs.fasta")
        # create a link btw. contigs.fasta -> Trinity.fasta
        relpth = os.path.relpath(trinity_fname, reads.r1.dir)
        os.symlink(relpth, contig_lname)
        # create a link btw. ../contigs/genus-species.contigs.fasta -> Trinity.fasta
        relpth = os.path.relpath(trinity_fname, contig_dir)
        contig_lname = os.path.join(contig_dir, sample)
        os.symlink(relpth, "{}.contigs.fasta".format(contig_lname))
    except:
        log.warn("Unable to symlink {} to {}".format(trinity_fname, contig_lname))


def cleanup_trinity_assembly_folder(pth, log):
    log.info("Removing extraneous Trinity files")
    files = glob.glob(os.path.join(pth, '*'))
    # check the names to make sure we're not deleting something improperly
    names = [os.path.basename(f) for f in files]
    try:
        assert "Trinity.fasta" in names
        assert "trinity.log" in names
    except:
        raise IOError("Neither Trinity.fasta nor trinity.log were found in output.")
    for file in files:
        if not os.path.basename(file) in ("Trinity.fasta", "trinity.log"):
            if os.path.isfile(file) or os.path.islink(file):
                os.remove(file)
            elif os.path.isdir(file):
                shutil.rmtree(file)


def main():
    # get args and options
    args = get_args()
    # setup logging
    log, my_name = setup_logging(args)
    # get the input data
    log.info("Getting input filenames and creating output directories")
    input = get_input_data(args.config, args.dir)
    # create the output directory if it does not exist
    if not os.path.isdir(args.output):
        os.makedirs(args.output)
    else:
        pass
    # make the symlink directory within the output directory
    contig_dir = os.path.join(args.output, 'contigs')
    if not os.path.isdir(contig_dir):
        os.makedirs(contig_dir)
    else:
        pass
    # Get path to trinity.  Standard name is `Trinity.pl`.
    # I usually symlink to `trinity`
    #TODO:  Change this to system "which" - this is just to flaky in certain cases
    try:
        trinity = which('trinity')[0]
    except EnvironmentError:
        trinity = which('Trinity.pl')[0]
    except:
        raise EnvironmentError("Cannot find Trinity.  Ensure it is installed and in your $PATH")
    for group in input:
        sample, dir = group
        # pretty print taxon status
        text = " Processing {} ".format(sample)
        log.info(text.center(65, "-"))
        # make a directory for sample-specific assemblies
        sample_dir = os.path.join(args.output, sample)
        os.makedirs(sample_dir)
        # determine how many files we're dealing with
        reads = get_input_files(dir, args.subfolder, log)
        # copy the read data over, combine singletons with read 1
        # and run the assembly for PE data.
        if reads.r1 and reads.r2 and reads.singleton:
            copy_read_data(reads, sample_dir, log)
            combine_read_data(reads, log)
            output = run_trinity_pe(trinity, reads, args.cores, log)
            if args.clean:
                cleanup_trinity_assembly_folder(output, log)
        # we don't need to combine singleton files here.  copy
        # the read data over and run the assembly for PE data
        elif reads.r1 and reads.r2:
            copy_read_data(reads, sample_dir, log)
            output = run_trinity_pe(trinity, reads, args.cores, log)
            if args.clean:
                cleanup_trinity_assembly_folder(output, log)
        # here, we don't have PE data, so copy the file over
        # and run the assembly for SE data
        elif reads.r1:
            copy_read_data(reads, sample_dir, log)
            output = run_trinity_se(trinity, reads, args.cores, log)
            if args.clean:
                cleanup_trinity_assembly_folder(output, log)
        # generate symlinks to assembled contigs
        generate_symlinks(contig_dir, sample, reads, log)
    text = " Completed {} ".format(my_name)
    log.info(text.center(65, "="))

if __name__ == '__main__':
    main()
