#!/usr/bin/env sh
# Copyright (C) 2019 The ESPResSo project
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

ISSUE_NUMBER=$(curl -s -G https://api.github.com/search/issues \
     --data-urlencode "q=\"CI failed for merged PR\" org:espressomd repo:espresso is:open is:issue in:title" \
     --data-urlencode "q=${CI_PIPELINE_ID} org:espressomd repo:espresso is:open is:issue in:body" | jq '.items[0] .number')

if [ "${ISSUE_NUMBER}" != "null" ]; then
    curl -s "https://api.github.com/repos/espressomd/espresso/issues/${ISSUE_NUMBER}" \
         -H "Accept: application/vnd.github.full+json" \
         -H "Content-Type: application/json" \
         -H "Authorization: token ${GITHUB_TOKEN}" \
         -X PATCH \
         -d "{\"state\": \"closed\" }"
fi
