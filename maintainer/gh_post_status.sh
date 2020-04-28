#!/usr/bin/env sh
# Copyright (C) 2017,2019 The ESPResSo project
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

[ "${#}" -eq 1 ] || exit -1

GIT_COMMIT=$(git rev-parse HEAD)
URL="https://gitlab.icp.uni-stuttgart.de/espressomd/espresso/pipelines/${CI_PIPELINE_ID}"
STATUS="${1}"
curl "https://api.github.com/repos/espressomd/espresso/statuses/${GIT_COMMIT}" \
     -H "Authorization: token ${GITHUB_TOKEN}" \
     -H "Content-Type: application/json" \
     -X POST \
     -d "{\"state\": \"${STATUS}\", \"context\": \"ICP GitLab CI\", \"target_url\": \"${URL}\"}"
