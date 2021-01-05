#
# Copyright (C) 2018-2020 The ESPResSo project
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

import os
import requests

if not os.environ.get('CI_COMMIT_REF_NAME', '').startswith('PR-'):
    raise RuntimeError("Not a pull request.")

PR = os.environ['CI_COMMIT_REF_NAME'][3:]
API = 'https://api.github.com'
URL = f'{API}/repos/espressomd/espresso/issues/{PR}/comments'
HEADERS = {'Authorization': 'token ' + os.environ['GITHUB_TOKEN']}
CI_JOB_URL = os.environ['CI_JOB_URL']
LOGIN = 'espresso-ci'


def delete_comments_by_token(token):
    """
    Delete posts made by the espresso-ci bot containing the string in ``token``.
    """
    comments = requests.get(URL, headers=HEADERS)
    comments.raise_for_status()
    for comment in comments.json():
        if comment['user']['login'] == LOGIN and token in comment['body']:
            print(f"deleting {comment['url']}")
            response = requests.delete(comment['url'], headers=HEADERS)
            response.raise_for_status()


def post_message(message):
    """
    Post a message using the espresso-ci bot.
    """
    response = requests.post(URL, headers=HEADERS, json={'body': message})
    response.raise_for_status()
