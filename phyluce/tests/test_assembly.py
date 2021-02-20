#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2018 Brant Faircloth || http://faircloth-lab.org/

All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 06 July 2018 15:27 CDT (-0500)
"""

import os
import re
import glob
import shutil
import platform
import subprocess

import pytest
from Bio import SeqIO

import pdb


@pytest.fixture(autouse=True)
def cleanup_files(request):
    """cleanup extraneous log files"""

    def clean():
        log_files = os.path.join(
            request.config.rootdir, "phyluce", "tests", "*.log"
        )
        for file in glob.glob(log_files):
            os.remove(file)

    request.addfinalizer(clean)


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


@pytest.fixture(scope="module")
def raw_dir(request):
    directory = os.path.join(
        request.config.rootdir,
        "phyluce",
        "tests",
        "test-expected",
        "raw-reads-short",
    )
    return directory


@pytest.fixture(scope="module")
def a_conf(request):
    return os.path.join(
        request.config.rootdir,
        "phyluce",
        "tests",
        "test-conf",
        "assembly-short.conf",
    )


def get_contig_lengths_and_counts(contigs):
    with open(contigs) as contig_file:
        contig_count = 0
        contig_length = 0
        for seq in SeqIO.parse(contig_file, "fasta"):
            contig_count += 1
            contig_length += len(seq)
    return contig_count, contig_length


def test_spades_assembly(o_dir, raw_dir, e_dir, request):
    program = "bin/assembly/phyluce_assembly_assemblo_spades"
    cmd = [
        os.path.join(request.config.rootdir, program),
        "--dir",
        raw_dir,
        "--cores",
        "1",
        "--output",
        "{}".format(os.path.join(o_dir, "spades")),
        "--log-path",
        o_dir,
    ]
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = proc.communicate()
    assert proc.returncode == 0, print("""{}""".format(stderr.decode("utf-8")))
    print(stderr)
    observed_count, observed_length = get_contig_lengths_and_counts(
        os.path.join(
            o_dir,
            "spades",
            "contigs",
            "alligator-mississippiensis.contigs.fasta",
        )
    )
    expected_count, expected_length = get_contig_lengths_and_counts(
        os.path.join(
            e_dir,
            "spades",
            "contigs",
            "alligator_mississippiensis.contigs.fasta",
        )
    )
    assert observed_count == expected_count
    assert observed_length == expected_length


def test_abyss_assembly(o_dir, raw_dir, e_dir, request):
    program = "bin/assembly/phyluce_assembly_assemblo_abyss"
    cmd = [
        os.path.join(request.config.rootdir, program),
        "--dir",
        raw_dir,
        "--cores",
        "1",
        "--output",
        "{}".format(os.path.join(o_dir, "abyss")),
        "--abyss-se",
        "--log-path",
        o_dir,
    ]
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = proc.communicate()
    assert proc.returncode == 0, print("""{}""".format(stderr.decode("utf-8")))
    observed_count, observed_length = get_contig_lengths_and_counts(
        os.path.join(
            o_dir,
            "abyss",
            "contigs",
            "alligator-mississippiensis.contigs.fasta",
        )
    )
    expected_count, expected_length = get_contig_lengths_and_counts(
        os.path.join(
            e_dir,
            "abyss",
            "contigs",
            "alligator_mississippiensis.contigs.fasta",
        )
    )
    assert observed_count == expected_count
    assert observed_length == expected_length


def test_velvet_assembly(o_dir, raw_dir, e_dir, request):
    program = "bin/assembly/phyluce_assembly_assemblo_velvet"
    cmd = [
        os.path.join(request.config.rootdir, program),
        "--dir",
        raw_dir,
        "--cores",
        "1",
        "--output",
        "{}".format(os.path.join(o_dir, "velvet")),
        "--log-path",
        o_dir,
    ]
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = proc.communicate()
    assert proc.returncode == 0, print("""{}""".format(stderr.decode("utf-8")))
    observed_count, observed_length = get_contig_lengths_and_counts(
        os.path.join(
            o_dir,
            "velvet",
            "contigs",
            "alligator-mississippiensis.contigs.fasta",
        )
    )
    expected_count, expected_length = get_contig_lengths_and_counts(
        os.path.join(
            e_dir,
            "velvet",
            "contigs",
            "alligator_mississippiensis.contigs.fasta",
        )
    )
    assert observed_count == expected_count
    # velvet results change all the time - lets make sure this is
    # simply almost equal
    assert abs(observed_length - expected_length) < 2000
