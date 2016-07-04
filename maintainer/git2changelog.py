#!/usr/bin/python
# Copyright (C) 2012,2013,2014,2015,2016 The ESPResSo project
# Copyright (C) 2011 Olaf Lenz
# Copyright 2008 Marcus D. Hanwell <marcus@cryos.org>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
from __future__ import print_function
import string, re, os

# Execute git log with the desired command line options.
fin = os.popen('git log --summary --stat --no-merges --date=short 3.0.1..', 'r')

# Set up the loop variables in order to locate the blocks we want
authorFound = False
dateFound = False
messageFound = False
filesFound = False
message = ""
messageNL = False
files = ""
prevAuthorLine = ""
commitId = ""

# The main part of the loop
for line in fin:
    # The commit line marks the start of a new commit object.
    m = re.match('^commit (.*)$', line)
    if m is not None:
        commitId = m.group(1)
        # Start all over again...
        authorFound = False
        dateFound = False
        messageFound = False
        messageNL = False
        message = ""
        filesFound = False
        files = ""
        continue

    # Match the author line and extract the part we want
    m = re.match('^Author:\s*(.*)\s*$', line)
    if m is not None:
        author = m.group(1)
        authorFound = True
        continue

    # Match the date line
    m = re.match('^Date:\s*(.*)\s*$', line)
    if m is not None:
        date = m.group(1)
        dateFound = True
        continue

    # The svn-id lines are ignored
    # The sign off line is ignored too
    if re.search('git-svn-id:|^Signed-off-by', line) >= 0:
        continue

    # Extract the actual commit message for this commit
    if not (authorFound & dateFound & messageFound):
        # Find the commit message if we can
        if len(line) == 1:
            if messageNL:
                messageFound = True
            else:
                messageNL = True
        elif len(line) == 4:
            messageFound = True
        else:
            if len(message) == 0:
                message = message + line.strip()
            else:
                message = message + " " + line.strip()

    # If this line is hit all of the files have been stored for this commit
    if re.search('files changed', line) >= 0:
        filesFound = True
        continue
    # Collect the files for this commit. FIXME: Still need to add +/- to files
    elif authorFound & dateFound & messageFound:
        fileList = re.split(' \| ', line, 2)
        if len(fileList) > 1:
            if len(files) > 0:
                files = files + ", " + fileList[0].strip()
            else:
                files = fileList[0].strip()

    # All of the parts of the commit have been found - write out the entry
    if authorFound & dateFound & messageFound & filesFound:
        # First the author line, only outputted if it is the first for that
        # author on this day
        authorLine = date + " " + author
        if len(prevAuthorLine) == 0:
            print(authorLine)
        elif authorLine == prevAuthorLine:
            pass
        else:
            print(("\n" + authorLine))

        # Assemble the actual commit message line(s) and limit the line length
        # to 80 characters.
        commitLine = "* " + files + ": " + message
        i = 0
        commit = ""
        while i < len(commitLine):
            if len(commitLine) < i + 78:
                commit = commit + "\n " + commitLine[i:len(commitLine)]
                break
            index = commitLine.rfind(' ', i, i+78)
            if index > i:
                commit = commit + "\n " + commitLine[i:index]
                i = index+1
            else:
                commit = commit + "\n " + commitLine[i:78]
                i = i+79

        # Write out the commit line
        print(commit)

        #Now reset all the variables ready for a new commit block.
        authorFound = False
        dateFound = False
        messageFound = False
        messageNL = False
        message = ""
        filesFound = False
        files = ""
        commitId = ""
        prevAuthorLine = authorLine

# Close the input and output lines now that we are finished.
fin.close()

