#!/usr/bin/env python
# encoding: utf-8
"""
File: assembly.py
Author: Brant Faircloth

Created by Brant Faircloth on 30 September 2013 10:09 PDT (-0700)
Copyright (c) 2013 Brant C. Faircloth. All rights reserved.

"""
import os
import re
import glob
import configparser

# import pdb


class Read(object):
    """Fastq reads"""

    def __init__(self, dir, file):
        self.dir = dir
        self.file = file
        if dir is not None and file is not None:
            self.pth = os.path.join(dir, file)
        else:
            self.pth = None

    def __str__(self):
        return "{} fastq read".format(self.file)

    def __repr__(self):
        return "<{}.{} instance at {}>".format(
            self.file, self.__class__.__name__, hex(id(self))
        )


class Fastqs(object):
    """Container for fastq data"""

    def __init__(self):
        self.r1 = None
        self.r2 = None
        self.singleton = None
        self.type = None
        self.gzip = False
        self.type = "fastq"
        self.reads = ()

    def __str__(self):
        return "Fastq container of R1, R2, Singletons"

    def set_read(self, read, dir, file):
        if read == "r1":
            self.r1 = Read(dir, file)
            self.reads += ((self.r1),)
        elif read == "r2":
            self.r2 = Read(dir, file)
            self.reads += ((self.r2),)
        elif read == "singleton":
            self.singleton = Read(dir, file)
            self.reads += ((self.singleton),)


class Fastas(Fastqs):
    """Container for fasta data"""

    def __init__(self):
        Fastqs.__init__(self)
        self.type = "fasta"


def check_for_fastq(dir, subfolder):
    types = (
        "*.fastq.gz",
        "*.fastq.gzip",
        "*.fq.gz",
        "*fq.gzip",
        "*.fq",
        "*.fastq",
    )
    files = []
    for type in types:
        files.extend(glob.glob(os.path.join(dir, subfolder, type)))
    return files


def check_for_fasta(dir, subfolder):
    types = (
        "*.fasta.gz",
        "*.fasta.gzip",
        "*.fa.gz",
        "*fa.gzip",
        "*.fa",
        "*.fasta",
    )
    files = []
    for type in types:
        files.extend(glob.glob(os.path.join(dir, subfolder, type)))
    return files


def get_input_files(dir, subfolder, log):
    log.info("Finding fastq/fasta files")
    fastq_files = check_for_fastq(dir, subfolder)
    fasta_files = check_for_fasta(dir, subfolder)
    if fastq_files and fasta_files:
        raise IOError("There are both fasta and fastq files in {}".format(dir))
    if not fastq_files and not fasta_files:
        raise IOError("There are not appropriate files in {}".format(dir))
    if fastq_files:
        log.info("File type is fastq")
        fq = Fastqs()
        files = fastq_files
    elif fasta_files:
        log.info("File type is fasta")
        fq = Fastas()
        files = fasta_files
    # get dirname of first file
    dir = os.path.dirname(files[0])
    ext = set()
    for f in files:
        # get file extension
        ext.add(os.path.splitext(f)[-1])
        # get file name
        fname = os.path.basename(f)
        # find which reach this is
        regex = (
            r"(?:.*)[_-](?:READ|Read|R)(\d)*[_-]*(singleton|unpaired)*(?:.*)"
        )
        match = re.search(
            regex,
            fname,
        )
        try:
            if match.groups()[0] == "1":
                assert fq.r1 is None
                fq.set_read("r1", dir, fname)
            elif match.groups()[0] == "2":
                assert fq.r2 is None
                fq.set_read("r2", dir, fname)
            elif (
                match.groups()[1] == "singleton"
                or match.groups()[1] == "unpaired"
            ):
                assert fq.singleton is None
                fq.set_read("singleton", dir, fname)
        except:
            raise IOError(
                "The appear to be multiple files for R1/R2/Singleton reads"
            )
    if len(ext) != 1:
        raise IOError("Files are of different types (e.g. gzip and fastq)")
    if ".gzip" in ext or ".gz" in ext:
        fq.gzip = True
    return fq


def get_input_data(config, dir):
    if config is not None:
        conf = configparser.ConfigParser()
        conf.optionxform = str
        conf.read(config)
        groups = conf.items("samples")
        # expand params when given userdirs or relative paths
        groups = [
            (g[0], os.path.abspath(os.path.expanduser(g[1]))) for g in groups
        ]
        for sample in groups:
            try:
                assert os.path.isdir(
                    os.path.abspath(os.path.expanduser(sample[1]))
                )
            except:
                raise IOError("{} is not a directory".format(sample[1]))
    else:
        groups = []
        for name in glob.glob(os.path.join(dir, "*")):
            groups.append((os.path.basename(name), name))
    return groups
