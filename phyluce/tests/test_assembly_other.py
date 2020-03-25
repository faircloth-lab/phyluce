#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2018 Brant Faircloth || http://faircloth-lab.org/

All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 06 July 2018 15:37 CDT (-0500)
"""

import os
import shutil
import subprocess
import configparser

# from phyluce.tests.common import get_contig_lengths_and_counts

import pytest

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


@pytest.fixture(scope="module")
def conf_dir(request):
    directory = os.path.join(
        request.config.rootdir, "phyluce", "tests", "test-conf"
    )
    return directory


def get_cmd(o_dir, e_dir, conf_dir, output_config, request, incomplete=False):
    program = "bin/assembly/phyluce_assembly_get_match_counts"
    cmd = [
        os.path.join(request.config.rootdir, program),
        "--locus-db",
        os.path.join(e_dir, "probe-match", "probe.matches.sqlite"),
        "--taxon-list-config",
        os.path.join(conf_dir, "taxon-set.conf"),
        "--taxon-group",
        "all",
        "--output",
        output_config,
        "--log-path",
        o_dir,
    ]
    if not incomplete:
        return cmd
    else:
        cmd.append("--incomplete-matrix")
        return cmd


def test_get_match_counts_complete(o_dir, e_dir, conf_dir, request):
    output_config = os.path.join(o_dir, "taxon-set.conf")
    cmd = get_cmd(o_dir, e_dir, conf_dir, output_config, request)
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = proc.communicate()
    obs_config = configparser.RawConfigParser(allow_no_value=True)
    obs_config.optionxform = str
    obs_config.read(output_config)
    expected_config = os.path.join(e_dir, "taxon-set.complete.conf")
    exp_config = configparser.RawConfigParser(allow_no_value=True)
    exp_config.optionxform = str
    exp_config.read(expected_config)
    assert obs_config == exp_config


def test_get_match_counts_incomplete(o_dir, e_dir, conf_dir, request):
    output_config = os.path.join(o_dir, "taxon-set.conf")
    cmd = get_cmd(
        o_dir, e_dir, conf_dir, output_config, request, incomplete=True
    )
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
    )
    stdout, stderr = proc.communicate()
    obs_config = configparser.RawConfigParser(allow_no_value=True)
    obs_config.optionxform = str
    obs_config.read(output_config)
    expected_config = os.path.join(e_dir, "taxon-set.incomplete.conf")
    exp_config = configparser.RawConfigParser(allow_no_value=True)
    exp_config.optionxform = str
    exp_config.read(expected_config)
    assert obs_config == exp_config


"""
def test_explode_get_fastas_file(o_dir):
    program = "bin/assembly/phyluce_assembly_explode_get_fastas_file"
    cmd = [
        os.path.join(ROOTDIR, program),
        "--config",
        a_conf,
        "--cores",
        "1",
        "--output",
        "{}".format(os.path.join(o_dir, "spades")),
        "--log-path",
        o_dir,
    ]
"""
