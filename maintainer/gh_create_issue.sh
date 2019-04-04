#!/usr/bin/env bash

URL=$(echo "https://gitlab.icp.uni-stuttgart.de/espressomd/espresso/pipelines/${CI_PIPELINE_ID}")

curl -s "https://api.github.com/repos/espressomd/espresso/issues" \
     -H "Accept: application/vnd.github.full+json" \
     -H "Content-Type: application/json" \
     -H "Authorization: token $GITHUB_TOKEN" \
     -X POST \
     -d "{\"title\": \"CI build failed for merged PR\", \"body\": \"$URL\" }"
