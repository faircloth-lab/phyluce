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
import re
import sys
import glob
import shutil
import logging
import argparse
import subprocess
import ConfigParser
from phyluce.helpers import FullPaths, is_dir, is_file
from phyluce.third_party import which

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


class Read():
    """Fastq reads"""
    def __init__(self, dir, file):
        self.dir = dir
        self.file = file

    def __str__(self):
        return "{} fastq read".format(self.read)

    def __repr__(self):
        return "<{}.{} instance at {}>".format(self.read, self.__class__.__name__, hex(id(self)))


class Fastqs():
    """Container for fastq data"""
    def __init__(self):
        self.r1 = None
        self.r2 = None
        self.singleton = None
        self.type = None
        self.reads = ()

    def __str__(self):
        return "Fastq container of R1, R2, Singletons"

    def set_read(self, read, dir, file):
        if read == 'r1':
            self.r1 = Read(dir, file)
            self.reads += ((self.r1),)
        elif read == 'r2':
            self.r2 = Read(dir, file)
            self.reads += ((self.r2),)
        elif read == 'singleton':
            self.singleton = Read(dir, file)
            self.reads += ((self.singleton),)


def get_fastq_input_files(dir, subfolder, log):
    log.info("Finding fastq files")
    types = ("*.fastq.gz", "*.fastq.gzip", "*.fq.gz", "*fq.gzip", "*.fq", "*.fastq")
    files = []
    for type in types:
        files.extend(glob.glob(os.path.join(dir, subfolder, type)))
    if not files:
        raise IOError("There are not appropriate files in {}".format(dir))
    fq = Fastqs()
    # get dirname of first file
    dir = os.path.dirname(files[0])
    ext = set()
    for f in files:
        # get file extension
        ext.add(os.path.splitext(f)[-1])
        # get file name
        fname = os.path.basename(f)
        # find which reach this is
        match = re.search("(?:.*)[_-](?:READ|Read|R)(\d)*[_-]*(singleton)*(?:.*)", fname)
        try:
            if match.groups()[0] == '1':
                assert fq.r1 is None
                fq.set_read('r1', dir, fname)
            elif match.groups()[0] == '2':
                assert fq.r2 is None
                fq.set_read('r2', dir, fname)
            elif match.groups()[1] == 'singleton':
                assert fq.singleton is None
                fq.set_read('singleton', dir, fname)
        except:
            raise IOError("The appear to be multiple files for R1/R2/Singleton reads")
    if len(ext) != 1:
        raise IOError("Files are of different types (e.g. gzip and fastq)")
    if '.gzip' in ext or '.gz' in ext:
        fq.type = 'gzip'
    else:
        fq.type = 'fastq'
    return fq


def get_input_data(config, dir):
    if config is not None:
        conf = ConfigParser.ConfigParser()
        conf.optionxform = str
        conf.read(config)
        groups = conf.items('samples')
        for sample in groups:
            try:
                assert os.path.isdir(sample[1])
            except:
                raise IOError("{} is not a directory".format(sample[1]))
    else:
        groups = []
        for name in glob.glob(os.path.join(dir, '*')):
            groups.append((os.path.basename(name), name))
    return groups


def copy_read_data(fastq, sample_dir, log):
    """We need to combine singleton reads with read1
    data to include those data so trinity can use them"""
    # copy read data to sample dir
    log.info("Copying raw read data to {}".format(sample_dir))
    for fq in fastq.reads:
        old_pth = os.path.join(fq.dir, fq.file)
        new_pth = os.path.join(sample_dir, fq.file)
        shutil.copy(old_pth, new_pth)
        # reset fq dirname
        fq.dir = sample_dir


def combine_read_data(fastq, log):
    log.info("Combining singleton reads with R1 data")
    # setup cat command for subprocess
    cmd = [
        "cat",
        os.path.join(fastq.singleton.dir, fastq.singleton.file)
    ]
    # cat the singleton fastq contents into the r1 fastq file
    with open(os.path.join(fastq.r1.dir, fastq.r1.file), 'a') as outf:
        proc = subprocess.Popen(cmd, stdout=outf)
        proc.communicate()
    # remove the singleton file
    os.remove(os.path.join(fastq.singleton.dir, fastq.singleton.file))
    # unset the singleton fastq record
    fastq.set_read('singleton', None, None)


def run_trinity_pe(trinity, fastq, cores, clean, log):
    log.info("Running Trinity.pl for PE data")
    cmd = [
        trinity,
        "--seqType",
        "fq",
        "--JM",
        "10G",
        "--left",
        os.path.join(fastq.r1.dir, fastq.r1.file),
        "--right",
        os.path.join(fastq.r2.dir, fastq.r2.file),
        "--CPU",
        str(cores),
        "--output",
        os.path.join(fastq.r1.dir)
    ]
    if clean:
        cmd.append("--full_cleanup")
    try:
        with open(os.path.join(fastq.r1.dir, 'trinity.log'), 'w') as outf:
            proc = subprocess.Popen(cmd, stdout=outf)
            proc.communicate()
    except:
        log.critical("Could not assemble {}".format(fastq.r1.dir))
    try:
        # trinity converts files to unzipped fastq. delete
        # gzips and fastas.  These data are also in `both.fa`
        # which is created at the start of Trinity.
        if fastq.type == 'gzip':
            r1 = os.path.join(fastq.r1.dir, fastq.r1.file)
            r2 = os.path.join(fastq.r2.dir, fastq.r2.file)
            r1_fastq = os.path.splitext(r1)[0]
            r2_fastq = os.path.splitext(r2)[0]
            for f in [r1, r2, r1_fastq, r2_fastq]:
                os.remove(f)
        else:
            r1 = os.path.join(fastq.r1.dir, fastq.r1.file)
            r2 = os.path.join(fastq.r2.dir, fastq.r2.file)
            for f in [r1, r2]:
                os.remove(f)
    except:
        log.warn("Did not clean all fastq files from {}".format(fastq.r1.dir))


def run_trinity_se(trinity, fastq, cores, clean, log):
    log.info("Running Trinity.pl for SE data")
    cmd = [
        trinity,
        "--seqType",
        "fq",
        "--JM",
        "10G",
        "--single",
        os.path.join(fastq.r1.dir, fastq.r1.file),
        "--CPU",
        str(cores),
        "--output",
        os.path.join(fastq.r1.dir)
    ]
    if clean:
        cmd.append("--full_cleanup")
    try:
        with open(os.path.join(fastq.r1.dir, 'trinity.log'), 'w') as outf:
            proc = subprocess.Popen(cmd, stdout=outf)
            proc.communicate()
    except:
        log.critical("Could not assemble {}".format(fastq.r1.dir))
    try:
        # trinity converts files to unzipped fastq. delete
        # gzips and fastas.  These data are also in `both.fa`
        # which is created at the start of Trinity.
        if fastq.type == 'gzip':
            r1 = os.path.join(fastq.r1.dir, fastq.r1.file)
            r1_fastq = os.path.splitext(r1)[0]
            for f in [r1, r1_fastq]:
                os.remove(f)
        else:
            r1 = os.path.join(fastq.r1.dir, fastq.r1.file)
            for f in [r1]:
                os.remove(f)
    except:
        log.warn("Did not clean all fastq files from {}".format(fastq.r1.dir))


def generate_symlinks(contig_dir, sample, fastq, clean, log):
    log.info("Symlinking assembled contigs into {}".format(contig_dir))
    pdb.set_trace()
    try:
        if not clean:
            trinity_fname = os.path.join(fastq.r1.dir, "Trinity.fasta")
            contig_lname = os.path.join(contig_dir, sample)
            os.symlink(trinity_fname, "{}.contigs.fasta".format(contig_lname))
        else:
            pass
    except:
        log.warn("Unable to symlink {} to {}".format(trinity_fname, contig_lname))


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
    log.info("=================== Starting ASSEMBLO-Trinity ===================")
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
        fastq = get_fastq_input_files(dir, args.subfolder, log)
        # copy the read data over, combine singletons with read 1
        # and run the assembly for PE data.
        if fastq.r1 and fastq.r2 and fastq.singleton:
            copy_read_data(fastq, sample_dir, log)
            combine_read_data(fastq, log)
            run_trinity_pe(trinity, fastq, args.cores, args.clean, log)
            #pdb.set_trace()
        # we don't need to combine singleton files here.  copy
        # the read data over and run the assembly for PE data
        elif fastq.r1 and fastq.r2:
            copy_read_data(fastq, sample_dir, log)
            run_trinity_pe(trinity, fastq, args.cores, args.clean, log)
        # here, we don't have PE data, so copy the file over
        # and run the assembly for SE data
        elif fastq.r1:
            copy_read_data(fastq, sample_dir, log)
            run_trinity_se(trinity, fastq, args.cores, args.clean, log)
        # generate symlinks to assembled contigs
        generate_symlinks(contig_dir, sample, fastq, args.clean, log)
    # pretty print some stuff
    #pdb.set_trace()


if __name__ == '__main__':
    main()
