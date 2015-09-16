#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
(c) 2015 Brant Faircloth || http://faircloth-lab.org/
All rights reserved.

This code is distributed under a 3-clause BSD license. Please see
LICENSE.txt for more information.

Created on 11 April 2015 16:57 CDT (-0500)
"""

import os
import sys
import ConfigParser

#import pdb

def get_user_path(program, binary, package_only=False):
    config = ConfigParser.ConfigParser()
    # make case sensitive
    config.optionxform = str
    if package_only:
        config.read(os.path.join(sys.prefix, 'config/phyluce.conf'))
    else:
        config.read([os.path.join(sys.prefix, 'config/phyluce.conf'), os.path.expanduser('~/.phyluce.conf')])
    # ensure program is in list
    pth = config.get(program, binary)
    # expand path as necessary - replace CONDA variable placeholder
    # with sys.prefix, otherwise default to normal path expansion
    if pth.startswith("$CONDA"):
        expand_pth = pth.replace("$CONDA", sys.prefix)
    else:
        expand_pth = os.path.abspath(os.path.expanduser(os.path.expandvars(pth)))
    return expand_pth


def get_user_param(section, param):
    config = ConfigParser.ConfigParser()
    # make case sensitive
    config.optionxform = str
    config.read([os.path.join(sys.prefix, 'config/phyluce.conf'), os.path.expanduser('~/.phyluce.conf')])
    return config.get(section, param)


def get_all_user_params(section):
    config = ConfigParser.ConfigParser()
    # make case sensitive
    config.optionxform = str
    config.read([os.path.join(sys.prefix, 'config/phyluce.conf'), os.path.expanduser('~/.phyluce.conf')])
    return [item[1] for item in config.items(section)]
