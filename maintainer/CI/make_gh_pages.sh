#!/usr/bin/env sh
#
# Copyright (C) 2021-2024 The ESPResSo project
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

set -e # exit on first non-zero return code
set -x # print each command line before executing it

if [ ! "${GITHUB_ACTIONS}" = "true" ]; then
  echo "This script is meant to be executed by a GitHub Actions workflow";
  exit 1;
fi

# checkout GitHub pages
cd "${HOME}"
git clone --quiet git@github.com:espressomd/espressomd.github.io.git
cd espressomd.github.io

# check if already up-to-date (i.e. the commit SHA of the
# generated docs is identical to the current commit SHA)
LAST_COMMIT=$(git log -1 --pretty='format:"%s"' remotes/origin/main)
NEXT_COMMIT="Documentation for ${GITHUB_SHA}"
if [ "${NEXT_COMMIT}" = "${LAST_COMMIT}" ]; then
  echo "Documentation already up-to-date.";
  exit 0;
fi

# download artifacts
gitlab_rest_api_get() {
  # shellcheck disable=SC2068
  curl --fail --location --request GET --header "PRIVATE-TOKEN: ${GITLAB_READ_API}" $@
}
rest_api="https://gitlab.icp.uni-stuttgart.de/api/v4"
project_id=390
# get the status of the most recent CI pipelines
gitlab_rest_api_get --silent "${rest_api}/projects/${project_id}/pipelines?per_page=60" > pipelines.json
# get the most recent successful nightly build (where tutorials are available)
# and the most recent successful build (where sphinx and doxygen are available)
pipeline_sched_id=$(jq '[ .[] | select(.ref=="python" and .status=="success" and .source=="schedule") ][0].id' pipelines.json)
pipeline_merge_id=$(jq '[ .[] | select(.ref=="python" and .status=="success") ][0].id' pipelines.json)
gitlab_rest_api_get --silent "${rest_api}/projects/${project_id}/pipelines/${pipeline_sched_id}/jobs?per_page=100" > jobs_scheduled_pipeline.json
gitlab_rest_api_get --silent "${rest_api}/projects/${project_id}/pipelines/${pipeline_merge_id}/jobs?per_page=100" > jobs_branch_head_pipeline.json
# get the tutorial, sphinx, and doxygen CI job ids and download their artifacts
tutorial_job_id=$(jq '.[] | select(.name=="run_tutorials").id' jobs_scheduled_pipeline.json)
doxygen_job_id=$(jq '.[] | select(.name=="run_doxygen").id' jobs_branch_head_pipeline.json)
sphinx_job_id=$(jq '.[] | select(.name=="check_sphinx").id' jobs_branch_head_pipeline.json)
gitlab_rest_api_get --output tutorials.zip "${rest_api}/projects/${project_id}/jobs/${tutorial_job_id}/artifacts"
gitlab_rest_api_get --output doxygen.zip   "${rest_api}/projects/${project_id}/jobs/${doxygen_job_id}/artifacts"
gitlab_rest_api_get --output sphinx.zip    "${rest_api}/projects/${project_id}/jobs/${sphinx_job_id}/artifacts"
if grep -F '<!DOCTYPE html>' ./*.zip; then
  echo "The artifacts could not be downloaded.";
  exit 1;
fi
unzip -q "*.zip"

# create a fresh main branch containing the docs of old releases
git config --global user.email "noreply@icp.uni-stuttgart.de"
git config --global user.name "espresso-ci"
git checkout -b new_main remotes/origin/releases

# generate the landing page by merging the branch containing the
# HTML theme and Markdown files, then convert them to HTML files
git merge --quiet --commit --no-edit --allow-unrelated-histories remotes/origin/landing_page
make tutorials.md
make videos.md
for filename in *.md; do
  make "${filename%.md}.html";
done
rm ./*_header.html
git add ./*.html
git clean -f
git rm ./*.md ./*.py Makefile
git commit --quiet -m "Generate landing page"

# add devel documentation
rsync -a --delete --exclude="*.md5" --exclude="*.map" "build/doc/doxygen/html/" dox/
rsync -a --delete --exclude=".buildinfo" "build/doc/sphinx/html/" doc/
rsync -a --delete "build/doc/tutorials/html/" tutorials/
git add doc dox tutorials
git commit --quiet -m "${NEXT_COMMIT}"

# deploy
git push -f origin new_main:main
