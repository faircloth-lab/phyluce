#!/usr/bin/env python
# encoding: utf-8

import os
import sys
import subprocess

__version__ = "1.7.2"
__default_config__ = os.path.join(sys.prefix, "phyluce/config/phyluce.conf")
__default_workflow_dir__ = os.path.join(sys.prefix, "phyluce/workflows")

# get a commit value, if possible.  if not running from git
# should default to static version
try:
    possible_git_dir = os.path.dirname(
        os.path.split(os.path.abspath(__file__))[0]
    )
    possible_git = os.path.join(possible_git_dir, ".git")
    cmd = ["git", "--git-dir", possible_git, "rev-parse", "--short", "HEAD"]
    proc = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
    )
    if proc.returncode == 0:
        stdout, stderr = proc.communicate()
        __git_commit__ = "{}".format(stdout.strip().decode("utf-8"))
    else:
        __git_commit__ = "None"

except FileNotFoundError:
    __git_commit__ = "None"

#
