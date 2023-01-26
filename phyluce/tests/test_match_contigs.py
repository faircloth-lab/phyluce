#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2020 Brant Faircloth || http://faircloth-lab.org/

All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 25 March 2020 10:24 CDT (-0500)
"""

import os
import re
import csv
import glob
import shutil
import platform
import subprocess

# from phyluce.tests.common import get_contig_lengths_and_counts

import pytest

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
def bait_pth(request):
    return os.path.join(
        request.config.rootdir,
        "phyluce",
        "tests",
        "probes",
        "uce-5k-probes.fasta",
    )


@pytest.fixture(scope="module")
def contig_pth(request, e_dir):
    return os.path.join(e_dir, "spades", "contigs")


def get_match_count_to_probes_results(pth):
    with open(pth, "r") as csvfile:
        reader = csv.reader(csvfile)
        # skip header
        next(reader)
        result = {
            item[0]: [item[1], item[2], item[3], item[4], item[5]]
            for item in reader
        }
    return result


def match_contigs_to_probes(o_dir, e_dir, bait_pth, contig_pth, request):
    out_pth = "{}".format(os.path.join(o_dir, "match_contigs"))
    csv_pth = "{}".format(
        os.path.join(o_dir, "match_contigs", "match_contig_results.csv")
    )
    program = "bin/assembly/phyluce_assembly_match_contigs_to_probes"
    cmd = [
        os.path.join(request.config.rootdir, program),
        "--probes",
        bait_pth,
        "--contigs",
        contig_pth,
        "--output",
        out_pth,
        "--log-path",
        o_dir,
        "--csv",
        csv_pth,
    ]
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = proc.communicate()
    assert proc.returncode == 0, print("""{}""".format(stderr.decode("utf-8")))
    print(stderr)
    observed_count = get_match_count_to_probes_results(csv_pth)
    expected_count = get_match_count_to_probes_results(
        os.path.join(e_dir, "probe-match", "probe_match_results.csv")
    )
    for k in observed_count.keys():
        assert observed_count[k] == expected_count[k]


def check_lastz_file(o_dir, e_dir):
    observed = os.path.join(
        o_dir, "match_contigs", "alligator_mississippiensis.contigs.lastz"
    )
    expected = os.path.join(
        e_dir, "probe-match", "alligator_mississippiensis.contigs.lastz"
    )
    with open(observed) as f1, open(expected) as f2:
        for ol, el in zip(f1, f2):
            assert ol == el


def test_match_contigs_to_probes(o_dir, e_dir, bait_pth, contig_pth, request):
    match_contigs_to_probes(o_dir, e_dir, bait_pth, contig_pth, request)
    check_lastz_file(o_dir, e_dir)


def test_match_contigs_to_barcodes(o_dir, e_dir, request):
    program = "bin/assembly/phyluce_assembly_match_contigs_to_barcodes"
    output = os.path.join(o_dir, "barcode-check")
    cmd = [
        os.path.join(request.config.rootdir, program),
        "--contigs",
        os.path.join(e_dir, "barcodes", "contigs"),
        "--barcodes",
        os.path.join(e_dir, "barcodes", "gallus.coi.fasta"),
        "--output",
        output,
    ]
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = proc.communicate()
    assert proc.returncode == 0, print("""{}""".format(stderr.decode("utf-8")))
    observed = (
        re.search("(Best.*\n)", stdout.decode("utf-8")).groups()[0].strip()
    )
    assert (
        "Best BOLD systems match for locus comp17283_c0_seq1: Anas poecilorhyncha [SIBJP030-10]"
        in observed
    )
