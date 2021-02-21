#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
(c) 2021 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 2021-02-18 T14:32:25-06:00
"""

import os
import csv
import glob
import shutil
import subprocess
import configparser

from phyluce.pth import get_user_path

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
def c_dir(request):
    directory = os.path.join(
        request.config.rootdir, "phyluce", "tests", "test-conf"
    )
    return directory


# this test always need to run first
def test_mapping_workflow(o_dir, e_dir, c_dir, request):
    program = "bin/workflow/phyluce_workflow"
    configfile = os.path.join(c_dir, "mapping.config.yaml")
    output = os.path.join(o_dir, "workflow-mapping")
    cmd = [
        os.path.join(request.config.rootdir, program),
        "--config",
        configfile,
        "--output",
        output,
        "--cores",
        "1",
        "--workflow",
        "mapping",
    ]
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = proc.communicate()
    assert proc.returncode == 0, print("""{}""".format(stderr.decode("utf-8")))
    observed_summary = os.path.join(
        output, "coverage", "all-taxon.summary.csv"
    )
    assert observed_summary, "Observed summary file missing"
    expected_summary = os.path.join(
        e_dir, "workflow-mapping.all-taxon.summary.csv"
    )
    observed = {}
    expected = {}
    with open(observed_summary) as observed_f:
        for line in observed_f:
            if line.startswith("taxon"):
                pass
            else:
                ls = line.strip().split(",")
                observed[ls[0]] = ls[1:]
    with open(expected_summary) as expected_f:
        for line in expected_f:
            if line.startswith("taxon"):
                pass
            else:
                ls = line.strip().split(",")
                expected[ls[0]] = ls[1:]
    for critter, values in observed.items():
        assert float(values[0]) == pytest.approx(
            float(expected[critter][0]), 1
        )
        assert float(values[1]) == pytest.approx(
            float(expected[critter][1]), 0.5
        )
        assert float(values[2]) == pytest.approx(
            float(expected[critter][2]), 0.5
        )


def test_correction_workflow(o_dir, e_dir, c_dir, request):
    program = "bin/workflow/phyluce_workflow"
    configfile = os.path.join(c_dir, "contig-correction.config.yaml")
    output = os.path.join(o_dir, "workflow-correction")
    cmd = [
        os.path.join(request.config.rootdir, program),
        "--config",
        configfile,
        "--output",
        output,
        "--cores",
        "1",
        "--workflow",
        "correction",
    ]
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = proc.communicate()
    assert proc.returncode == 0, print("""{}""".format(stderr.decode("utf-8")))
    output_files = glob.glob(os.path.join(output, "consensus", "*"))
    assert output_files, "There are no output files"
    for output_file in output_files:
        name = os.path.basename(output_file)
        expected_file = os.path.join(e_dir, "workflow-contig-correction", name)
        observed = SeqIO.to_dict(SeqIO.parse(output_file, "fasta"))
        expected = SeqIO.to_dict(SeqIO.parse(expected_file, "fasta"))
        for name, observed in observed.items():
            assert expected[name].seq == observed.seq


def test_correction_phasing(o_dir, e_dir, c_dir, request):
    program = "bin/workflow/phyluce_workflow"
    configfile = os.path.join(c_dir, "phasing.config.yaml")
    output = os.path.join(o_dir, "workflow-phasing")
    cmd = [
        os.path.join(request.config.rootdir, program),
        "--config",
        configfile,
        "--output",
        output,
        "--cores",
        "1",
        "--workflow",
        "phasing",
    ]
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = proc.communicate()
    assert proc.returncode == 0, print("""{}""".format(stderr.decode("utf-8")))
    output_files = glob.glob(os.path.join(output, "fastas", "*.fasta"))
    assert output_files, "There are no output files"
    # The comparison below can be inconsistent across platforms
    # because the 0 and the 1 haplotypes get shuffled around btw.
    # the expected and observed data.
    """
    for output_file in output_files:
        name = os.path.basename(output_file)
        expected_file = os.path.join(e_dir, "workflow-phasing", name)
        observed = SeqIO.to_dict(SeqIO.parse(output_file, "fasta"))
        expected = SeqIO.to_dict(SeqIO.parse(expected_file, "fasta"))
        for name, observed in observed.items():
            assert str(expected[name].seq) == str(observed.seq)
    """
