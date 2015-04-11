#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2014 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 26 June 2014 17:16 PDT (-0700)
"""

import os
import subprocess

from phyluce.pth import get_user_path


def index(log, sample, sample_dir, bam):
    log.info("Indexing BAM for {}".format(sample))
    cmd = [
        get_user_path("samtools", "samtools"),
        "index",
        bam
    ]
    samtools_out_fname = os.path.join(sample_dir, '{}.samtools-index-out.log'.format(sample))
    with open(samtools_out_fname, 'w') as samtools_out:
        proc = subprocess.Popen(cmd, stdout=samtools_out, stderr=subprocess.STDOUT)
        proc.communicate()


def create_faidx(log, sample, sample_dir, fasta):
    log.info("Indexing fasta for {}".format(sample))
    cmd = [
        get_user_path("samtools", "samtools"),
        "faidx",
        fasta
    ]
    samtools_out_fname = os.path.join(sample_dir, '{}.samtools-faidx-out.log'.format(sample))
    with open(samtools_out_fname, 'w') as samtools_out:
        proc = subprocess.Popen(cmd, stdout=samtools_out, stderr=subprocess.STDOUT)
        proc.communicate()
