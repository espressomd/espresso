#!/usr/bin/python
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

from __future__ import print_function
import os
import requests
import subprocess

if not os.environ['CI_COMMIT_REF_NAME'].startswith('PR-'):
    exit(0)

PR = os.environ['CI_COMMIT_REF_NAME'][3:]
URL = 'https://api.github.com/repos/espressomd/espresso/issues/' + \
      PR + '/comments?access_token=' + os.environ['GITHUB_TOKEN']
SIZELIMIT = 10000

# Delete all existing comments
comments = requests.get(URL)
for comment in comments.json():
    if comment['user']['login'] == 'espresso-ci' and 'style.patch' in comment['body']:
        requests.delete(comment['url'] + '?access_token=' +
                        os.environ['GITHUB_TOKEN'])

# If the working directory is not clean, post a new comment
if subprocess.call(["git", "diff-index", "--quiet", "HEAD", "--"]) != 0:
    comment = 'Your pull request does not meet our code formatting rules. '
    patch = subprocess.check_output(['git', '--no-pager', 'diff'])
    if len(patch) <= SIZELIMIT:
        comment += 'Specifically, I suggest you make the following changes:\n'
        comment += '```diff\n'
        comment += patch.replace('`', r'\`').strip()
        comment += '\n```\n'
        comment += 'To apply these changes, please do one of the following:\n'
    else:
        comment += 'To fix this, please do one of the following:\n'
    comment += '- You can download a patch with my suggested changes '
    comment += '[here](' + os.environ['CI_JOB_URL'] + \
               '/artifacts/raw/style.patch), '
    comment += 'inspect it and make changes manually.\n'
    comment += '- You can directly apply it to your repository by running '
    comment += '`curl ' + os.environ['CI_JOB_URL'] + \
               '/artifacts/raw/style.patch | git apply -`.\n'
    comment += '- You can run `maintainer/CI/fix_style.sh` to automatically fix your coding style. This is the same command that I have executed to generate the patch above, but it requires certain tools to be installed on your computer.\n\n'
    comment += 'You can run `gitlab-runner exec docker style` afterwards to check if your changes worked out properly.\n\n'
    comment += 'Please note that there are often multiple ways to correctly format code. As I am just a robot, I sometimes fail to identify the most aesthetically pleasing way. So please look over my suggested changes and adapt them where the style does not make sense.'

    if len(patch) > 0:
        requests.post(URL, json={'body': comment})
