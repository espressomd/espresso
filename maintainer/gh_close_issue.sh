#!/usr/bin/env bash

ISSUE_NUMBER=$(cat ${CI_PROJECT_DIR}/ISSUE_NUMBER.txt)
curl -s "https://api.github.com/repos/espressomd/espresso/issues/$ISSUE_NUMBER" \
     -H "Accept: application/vnd.github.full+json" \
     -H "Content-Type: application/json" \
     -H "Authorization: token $GITHUB_TOKEN" \
     -X PATCH \
     -d "{\"state\": \"closed\" }"
