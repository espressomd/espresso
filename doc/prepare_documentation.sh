#!/bin/sh

set -e;


if [ "${TRAVIS_PULL_REQUEST}" != "false" ]; then
    echo "INFO: This is a PR.";
    echo "INFO: Not building docs.";
    exit 0;
fi;


if [ "${TRAVIS_BRANCH}" != "master" ]; then
    echo "INFO: We are not on the master branch.";
    echo "INFO: Not building docs.";
    exit 0;
fi;


if [ -z "$GH_TOKEN" ]; then
    echo "INFO: The GitHub access token is not set.";
    echo "INFO: Not building docs.";
    exit 0;
fi;


if [ -z "$(git ls-remote --heads https://github.com/${TRAVIS_REPO_SLUG} gh-pages)" ]; then
    echo "INFO: The branch gh-pages does not exist.";
    echo "INFO: Not building docs.";
    exit 0;
fi;


# Set up git in documentation directory
mkdir -p "doc/";
cd "doc/";

git init;
git remote add origin https://github.com/${TRAVIS_REPO_SLUG};
git checkout -b gh-pages;

git config user.name "Travis CI";
git config user.email "travis@travis-ci.org";
