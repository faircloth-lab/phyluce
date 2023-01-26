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

from phyluce import __version__, __git_commit__

# import pdb


class ColorFormatter(logging.Formatter):
    """Logging Formatter to add colors and count warning / errors.

    Taken from Stackeverflow user Sergey Pleshakov from this URL
    https://stackoverflow.com/a/56944256"""

    grey = "\x1b[38;21m"
    yellow = "\x1b[33;21m"
    red = "\x1b[31;21m"
    bold_red = "\x1b[31;1m"
    reset = "\x1b[0m"
    format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"

    FORMATS = {
        logging.DEBUG: grey + format + reset,
        logging.INFO: grey + format + reset,
        logging.WARNING: yellow + format + reset,
        logging.ERROR: red + format + reset,
        logging.CRITICAL: bold_red + format + reset,
    }

    def format(self, record):
        log_fmt = self.FORMATS.get(record.levelno)
        formatter = logging.Formatter(log_fmt)
        return formatter.format(record)


def setup_logging(args):
    import __main__ as main

    my_name = os.path.basename(os.path.splitext(main.__file__)[0])
    log = logging.getLogger(my_name)
    console = logging.StreamHandler(sys.stdout)
    if args.log_path is not None:
        logfile = logging.FileHandler(
            os.path.join(args.log_path, "{}.log".format(my_name)),
            encoding="utf8",
        )
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
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    console.setFormatter(ColorFormatter())
    logfile.setFormatter(formatter)
    log.addHandler(console)
    log.addHandler(logfile)
    text = " Starting {} ".format(my_name)
    log.info(text.center(65, "="))
    log.info("Version: {}".format(str(__version__)))
    log.info("Commit: {}".format(str(__git_commit__)))
    for arg, value in sorted(vars(args).items()):
        log.info("Argument --{}: {}".format(arg, value))
    return log, my_name
