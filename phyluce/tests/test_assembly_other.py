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

# import subprocess

# from phyluce.tests.common import get_contig_lengths_and_counts

import pytest

# import pdb


@pytest.fixture(scope="module")
def o_dir(request):
    directory = os.path.join(request.config.rootdir, "test")
    os.mkdir(directory)

    def clean():
        shutil.rmtree(directory)

    request.addfinalizer(clean)
    return directory


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
