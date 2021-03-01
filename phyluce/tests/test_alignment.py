#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
(c) 2021 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 2021-02-13 T16:12:15-06:00
"""

import os
import re
import glob
import shutil
import subprocess

import pytest
from Bio import AlignIO

import pdb


@pytest.fixture(scope="module")
def o_dir(request):
    directory = os.path.join(
        request.config.rootdir, "phyluce", "tests", "test-observed"
    )
    os.mkdir(directory)

    def clean():
        shutil.rmtree(directory)

    request.addfinalizer(clean)
    return directory


@pytest.fixture(scope="module")
def e_dir(request):
    directory = os.path.join(
        request.config.rootdir, "phyluce", "tests", "test-expected"
    )
    return directory


def test_seqcap_align_mafft_untrim(o_dir, e_dir, request):
    program = "bin/align/phyluce_align_seqcap_align"
    output = os.path.join(o_dir, "mafft")
    cmd = [
        os.path.join(request.config.rootdir, program),
        "--input",
        os.path.join(e_dir, "taxon-set.incomplete.fasta"),
        "--output",
        output,
        "--taxa",
        "4",
        "--aligner",
        "mafft",
        "--output-format",
        "nexus",
        "--no-trim",
        "--cores",
        "1",
    ]
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = proc.communicate()
    assert proc.returncode == 0, print("""{}""".format(stderr.decode("utf-8")))
    output_files = glob.glob(os.path.join(output, "*"))
    assert output_files, "There are no output files"
    for output_file in output_files:
        name = os.path.basename(output_file)
        expected_file = os.path.join(e_dir, "mafft-no-trim", name)
        observed = open(output_file).read()
        expected = open(expected_file).read()
        assert observed == expected


def test_seqcap_align_muscle_untrim(o_dir, e_dir, request):
    program = "bin/align/phyluce_align_seqcap_align"
    output = os.path.join(o_dir, "muscle")
    cmd = [
        os.path.join(request.config.rootdir, program),
        "--input",
        os.path.join(e_dir, "taxon-set.incomplete.fasta"),
        "--output",
        output,
        "--taxa",
        "4",
        "--aligner",
        "muscle",
        "--output-format",
        "nexus",
        "--no-trim",
        "--cores",
        "1",
    ]
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = proc.communicate()
    assert proc.returncode == 0, print("""{}""".format(stderr.decode("utf-8")))
    output_files = glob.glob(os.path.join(output, "*"))
    assert output_files, "There are no output files"
    for output_file in output_files:
        name = os.path.basename(output_file)
        expected_file = os.path.join(e_dir, "muscle-no-trim", name)
        observed = open(output_file).read()
        expected = open(expected_file).read()
        assert observed == expected
