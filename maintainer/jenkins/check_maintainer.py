#!/usr/bin/python
from __future__ import print_function
import sys, os, subprocess

# Test for precious files in commit
current_commit = os.environ["GIT_COMMIT"]

if "GIT_PREVIOUS_COMMIT" in os.environ:
    prev_commit = os.environ["GIT_PREVIOUS_COMMIT"]
    refspec = "{}..{}".format(current_commit, prev_commit)
else:
    refspec = current_commit

gitout = subprocess.check_output(["git", "diff", "--name-only", refspec])
files = gitout.split()

warnings = []
for fname in files:
    if fname.startswith('maintainer/') or \
          fname == "src/myconfig-default.hpp":
    print("check_maintainer.py: Precious file was modified: {}".format(fname))
