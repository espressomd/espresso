#!/usr/bin/env bash

ISSUE_NUMBER=$(curl -s -G https://api.github.com/search/issues \
     --data-urlencode "q=\"CI failed for merged PR\" org:espressomd repo:espresso is:open is:issue in:title" \
     --data-urlencode "q=${CI_PIPELINE_ID} org:espressomd repo:espresso is:open is:issue in:body" | jq '.items[0] .number')

if [ "$ISSUE_NUMBER" != "null" ]; then
curl -s "https://api.github.com/repos/espressomd/espresso/issues/$ISSUE_NUMBER" \
     -H "Accept: application/vnd.github.full+json" \
     -H "Content-Type: application/json" \
     -H "Authorization: token $GITHUB_TOKEN" \
     -X PATCH \
     -d "{\"state\": \"closed\" }"
fi
