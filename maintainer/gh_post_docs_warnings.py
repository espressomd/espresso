#!/usr/bin/env python3
#
# Copyright (C) 2018-2019 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

import re
import os
import sys
import requests
import subprocess

if not os.environ['CI_COMMIT_REF_NAME'].startswith('PR-'):
    exit(0)

PR = os.environ['CI_COMMIT_REF_NAME'][3:]
URL = 'https://api.github.com/repos/espressomd/espresso/issues/' + \
      PR + '/comments?access_token=' + os.environ['GITHUB_TOKEN']
SIZELIMIT = 5000

doc_type, has_warnings, filepath_warnings = sys.argv[-3:]
has_warnings = has_warnings != '0'
prefix = {'sphinx': 'doc', 'doxygen': 'dox'}[doc_type]

# Delete all existing comments
comments = requests.get(URL)
for comment in comments.json():
    if comment['user']['login'] == 'espresso-ci' and prefix + '_warnings.sh' in comment['body']:
        requests.delete(comment['url'] + '?access_token=' +
                        os.environ['GITHUB_TOKEN'])

# If documentation raised warnings, post a new comment
if has_warnings:
    with open(filepath_warnings) as f:
        warnings = f.read().strip()
    warnings = warnings.replace('@', '\\')
    assert warnings.count('\n') >= 1, 'list of warnings is missing newlines'
    # the logfile must be guarded by backticks
    backticks = max(['``'] + re.findall('`+', warnings), key=len) + '`'
    assert len(backticks) < 12, 'cannot guard logfile warnings with backticks'
    # format message
    summary, warnings = warnings.split('\n', 1)
    comment = 'Your pull request does not meet our code documentation rules. '
    comment += summary + '\n\n' + backticks + '\n'
    if len(warnings) > SIZELIMIT:
        for line in warnings.split('\n'):
            if len(comment) + len(line) > SIZELIMIT - 200:
                break
            comment += line + '\n'
        comment = comment.rstrip() + '\n' + backticks + '\n'
        comment += (
            '\nThis list was truncated, check the [container logfile]'
            '({}) for the complete list.\n'.format(os.environ['CI_JOB_URL']))
    else:
        comment += warnings.rstrip() + '\n' + backticks + '\n'
    comment += (
        '\nYou can generate these warnings with `make -t; make {}; '
        '../maintainer/CI/{}_warnings.sh` using the maxset config. This is '
        'the same command that I have executed to generate the log above.'
        .format(doc_type, prefix))

    requests.post(URL, json={'body': comment})
