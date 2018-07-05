#!/usr/bin/env python
# encoding: utf-8

import os
import subprocess

__version__ = "1.6.4"

# get a dynamic version number, if possible.  if not running from git
# should default to static version
try:
    possible_git_dir = os.path.dirname(os.path.split(os.path.abspath(__file__))[0])
    possible_git = os.path.join(possible_git_dir, ".git")
    cmd = [
        "git",
        "--git-dir",
        possible_git,
        "rev-parse",
        "--short",
        "HEAD"
    ]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout, stderr = proc.communicate()
    __version__ = "git {}".format(stdout.strip())

except:
    pass
