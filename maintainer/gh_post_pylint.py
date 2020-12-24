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

if not os.environ.get('CI_COMMIT_REF_NAME', '').startswith('PR-'):
    print("Not a pull request. Exiting now.")
    exit(0)

import re
import sys
import gh_post

SIZELIMIT = 5000
TOKEN_ESPRESSO_CI = 'Pylint summary'

n_warnings, filepath_warnings = sys.argv[-2:]

# Delete obsolete posts
gh_post.delete_comments_by_token(TOKEN_ESPRESSO_CI)

# If pylint raised errors, post a new comment
if n_warnings != '0':
    with open(filepath_warnings) as f:
        warnings = f.read().strip()
    # the logfile must be guarded by backticks
    backticks = max(['``'] + re.findall('`+', warnings), key=len) + '`'
    assert len(backticks) < 12, 'cannot guard logfile warnings with backticks'
    # format message
    comment = 'Your pull request does not meet our code style rules. '
    comment += 'Pylint summary:\n' + backticks + '\n'
    if len(warnings) > SIZELIMIT:
        for line in warnings.split('\n'):
            if len(comment) + len(line) > SIZELIMIT - 200:
                break
            comment += line + '\n'
        comment = comment.rstrip() + '\n' + backticks + '\n'
        comment += (
            f'\nThis list was truncated, check the [container logfile]'
            f'({gh_post.CI_JOB_URL}) for the complete list.\n')
    else:
        comment += warnings.rstrip() + '\n' + backticks + '\n'
    comment += (
        '\nYou can generate these warnings with `maintainer/CI/fix_style.sh`. '
        'This is the same command that I have executed to generate the log above.'
    )
    assert TOKEN_ESPRESSO_CI in comment
    gh_post.post_message(comment)
