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

import pdb

def setup_logging(args):
    import __main__ as main
    import __init__ as init
    my_name = os.path.basename(os.path.splitext(main.__file__)[0])
    log = logging.getLogger(my_name)
    console = logging.StreamHandler(sys.stdout)
    if args.log_path is not None:
        logfile = logging.FileHandler(os.path.join(args.log_path, "{}.log".format(my_name)))
    else:
        logfile = logging.FileHandler("{}.log".format(my_name))
    if args.verbosity == "INFO":
        log.setLevel(logging.INFO)
        console.setLevel(logging.INFO)
        logfile.setLevel(logging.INFO)
    if args.verbosity == "WARN":
        log.setLevel(logging.WARN)
        console.setLevel(logging.WARN)
        logfile.setLevel(logging.WARN)
    if args.verbosity == "CRITICAL":
        log.setLevel(logging.CRITICAL)
        console.setLevel(logging.CRITICAL)
        logfile.setLevel(logging.CRITICAL)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    console.setFormatter(formatter)
    logfile.setFormatter(formatter)
    log.addHandler(console)
    log.addHandler(logfile)
    text = " Starting {} ".format(my_name)
    log.info(text.center(65, "="))
    log.info("Version: {}".format(init.__version__))
    for arg, value in sorted(vars(args).items()):
        log.info("Argument --{}: {}".format(arg, value))
    return log, my_name
