#!/usr/bin/env sh
#
# Copyright (C) 2019-2022 The ESPResSo project
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

URL="https://gitlab.icp.uni-stuttgart.de/espressomd/espresso/-/pipelines/${CI_PIPELINE_ID}"

curl -s "https://api.github.com/repos/espressomd/espresso/issues" \
     -H "Accept: application/vnd.github.full+json" \
     -H "Content-Type: application/json" \
     -H "Authorization: token ${GITHUB_TOKEN}" \
     -X POST \
     -d "{\"title\": \"CI build failed for merged PR\", \"body\": \"${URL}\" }"
