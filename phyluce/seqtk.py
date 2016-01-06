#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2015 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 28 December 2015 16:28 CST (-0600)
"""


import subprocess

from phyluce.pth import get_user_path

#import pdb


def fq_to_fa(log, sample, sample_dir, fastq, phase=None):
    if phase is None:
        log.info("Converting --Unphased-- FASTQ files to FASTA files")
    else:
        log.info("Creating REF/ALT allele FASTA file {0} from FASTQ {0}".format(phase))
    cmd = [
        get_user_path("seqtk", "seqtk"),
        "seq",
        "-a",
        fastq
    ]
    seqtk_out_fname = "{}.seqtk-seq-out.log".format(sample_dir)
    if phase is None:
        seqtk_fasta_fname = "{}.fasta".format(sample_dir)
    else:
        seqtk_fasta_fname = "{}.{}.fasta".format(sample_dir, phase)
    with open(seqtk_out_fname, 'w') as seqtk_out:
        with open(seqtk_fasta_fname, 'w') as seqtk_fasta:
            proc = subprocess.Popen(cmd, stdout=seqtk_fasta, stderr=seqtk_out)
            proc.communicate()
    return seqtk_fasta_fname
