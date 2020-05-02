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

import os
import requests
import subprocess

if not os.environ['CI_COMMIT_REF_NAME'].startswith('PR-'):
    exit(0)

PR = os.environ['CI_COMMIT_REF_NAME'][3:]
URL = 'https://api.github.com/repos/espressomd/espresso/issues/' + \
      PR + '/comments'
HEADERS = {'Authorization': 'token ' + os.environ['GITHUB_TOKEN']}
SIZELIMIT = 10000
TOKEN_ESPRESSO_CI = 'style.patch'

# Delete all existing comments
comments = requests.get(URL, headers=HEADERS)
comments.raise_for_status()
for comment in comments.json():
    if comment['user']['login'] == 'espresso-ci' and \
            TOKEN_ESPRESSO_CI in comment['body']:
        response = requests.delete(comment['url'], headers=HEADERS)
        response.raise_for_status()

MESSAGE = '''Your pull request does not meet our code formatting \
rules. {header}, please do one of the following:

- You can download a patch with my suggested changes \
  [here]({url}/artifacts/raw/style.patch), inspect it and make \
  changes manually.
- You can directly apply it to your repository by running \
  `curl {url}/artifacts/raw/style.patch | git apply -`.
- You can run `maintainer/CI/fix_style.sh` to automatically fix your coding \
  style. This is the same command that I have executed to generate the patch \
  above, but it requires certain tools to be installed on your computer.

You can run `gitlab-runner exec docker style` afterwards to check if your \
changes worked out properly.

Please note that there are often multiple ways to correctly format code. \
As I am just a robot, I sometimes fail to identify the most aesthetically \
pleasing way. So please look over my suggested changes and adapt them \
where the style does not make sense.\
'''

# If the working directory is not clean, post a new comment
if subprocess.call(["git", "diff-index", "--quiet", "HEAD", "--"]) != 0:
    patch = subprocess.check_output(['git', '--no-pager', 'diff'])
    if len(patch) <= SIZELIMIT:
        comment = 'Specifically, I suggest you make the following changes:'
        comment += '\n```diff\n'
        comment += patch.decode('utf-8').replace('`', r'\`').strip()
        comment += '\n```\n'
        comment += 'To apply these changes'
    else:
        comment = 'To fix this'
    message = MESSAGE.format(header=comment, url=os.environ['CI_JOB_URL'])

    if patch:
        assert TOKEN_ESPRESSO_CI in message
        response = requests.post(URL, headers=HEADERS,
                                 json={'body': message})
        response.raise_for_status()
