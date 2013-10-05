#!/usr/bin/env python
# encoding: utf-8
"""
File: logging.py
Author: Brant Faircloth

Created by Brant Faircloth on 03 October 2013 14:10 PDT (-0700)
Copyright (c) 2013 Brant C. Faircloth. All rights reserved.

Description:

"""

import os
import sys
import logging

#import pdb

def setup_logging(level, pth=None):
    import __main__ as main
    my_name = os.path.basename(os.path.splitext(main.__file__)[0])
    log = logging.getLogger(my_name)
    console = logging.StreamHandler(sys.stdout)
    if pth is not None:
        logfile = logging.FileHandler(os.path.join(pth, "{}.log".format(my_name)))
    else:
        logfile = logging.FileHandler("{}.log".format(my_name))
    if level == "INFO":
        log.setLevel(logging.INFO)
        console.setLevel(logging.INFO)
        logfile.setLevel(logging.INFO)
    if level == "WARN":
        log.setLevel(logging.WARN)
        console.setLevel(logging.WARN)
        logfile.setLevel(logging.WARN)
    if level == "CRITICAL":
        log.setLevel(logging.CRITICAL)
        console.setLevel(logging.CRITICAL)
        logfile.setLevel(logging.CRITICAL)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    console.setFormatter(formatter)
    logfile.setFormatter(formatter)
    log.addHandler(console)
    log.addHandler(logfile)
    return log, my_name
