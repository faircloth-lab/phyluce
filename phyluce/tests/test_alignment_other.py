#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
(c) 2021 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 2021-02-14 T10:44:15-06:00
"""


import os
import re
import glob
import shutil
import platform
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

    # def clean():
    #    shutil.rmtree(directory)

    # request.addfinalizer(clean)
    return directory


@pytest.fixture(scope="module")
def e_dir(request):
    directory = os.path.join(
        request.config.rootdir, "phyluce", "tests", "test-expected"
    )
    return directory


@pytest.mark.skipif(
    platform.processor() == "arm64", reason="Wont run on arm64"
)
def test_align_gblocks_trim(o_dir, e_dir, request):
    program = (
        "bin/align/phyluce_align_get_gblocks_trimmed_alignments_from_untrimmed"
    )
    output = os.path.join(o_dir, "mafft-gblocks")
    cmd = [
        os.path.join(request.config.rootdir, program),
        "--alignments",
        os.path.join(e_dir, "mafft"),
        "--output",
        output,
        "--input-format",
        "fasta",
        "--output-format",
        "nexus",
        "--cores",
        "1",
    ]
    pdb.set_trace()
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = proc.communicate()
    for output_file in glob.glob(os.path.join(output, "*")):
        name = os.path.basename(output_file)
        expected_file = os.path.join(e_dir, "mafft-gblocks", name)
        observed = open(output_file).read()
        expected = open(expected_file).read()
        assert observed == expected
