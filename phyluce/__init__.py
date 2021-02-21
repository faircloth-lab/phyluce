#!/usr/bin/env python
# encoding: utf-8

import os
import sys
import subprocess

static_version = "2.0"

__default_config__ = os.path.join(sys.prefix, "config/phyluce.conf")
__default_workflow_dir__ = os.path.join(sys.prefix, "workflows")

# get a dynamic version number, if possible.  if not running from git
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
    stdout, stderr = proc.communicate()
    __version__ = "{}_r_{}".format(
        static_version, stdout.strip().decode("utf-8")
    )
except:
    __version__ = static_version
