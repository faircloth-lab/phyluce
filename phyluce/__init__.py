import os
import subprocess

# get a dynamic version number, if possible.  if not running from git
# should default to static version
cwd = os.getcwd()
try:
    location = os.path.split(os.path.abspath(__file__))[0]
    os.chdir(location)
    cmd = [
        "git",
        "rev-parse",
        "--short",
        "HEAD"
    ]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    stdout, stderr = proc.communicate()
    if stdout.startswith("fatal:"):
        raise IOError("{}".format(stdout.strip()))
    else:
        __version__ = "git {}".format(stdout.strip())
    os.chdir(cwd)
except OSError:
    __version__ = "1.5.0"
    if not os.getcwd == cwd:
        os.chdir(cwd)
except IOError:
    __version__ = "1.5.0"
    if not os.getcwd == cwd:
        os.chdir(cwd)

