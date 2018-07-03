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


if [ -z "${GH_TOKEN}" ]; then
    echo "INFO: The GitHub access token is not set.";
    echo "INFO: Not building docs.";
    exit 0;
fi;


if [ -z "$(git ls-remote --heads https://github.com/${TRAVIS_REPO_SLUG} gh-pages)" ]; then
    echo "INFO: The branch gh-pages does not exist.";
    echo "INFO: Not building docs.";
    exit 0;
fi;


function gh_pages_prepare()
{
    mkdir -p "${TRAVIS_BUILD_DIR}/build/doc/html";
    cd "${TRAVIS_BUILD_DIR}/build/doc/html";

    git init;
    git remote add origin https://github.com/${TRAVIS_REPO_SLUG};
    git checkout -b gh-pages;

    git config user.name "Travis CI";
    git config user.email "travis@travis-ci.org";
}


function gh_pages_generate()
{
    cd "${TRAVIS_BUILD_DIR}/build";

    cmake ..;

    make doxygen;
}


function gh_pages_update()
{
    cd "${TRAVIS_BUILD_DIR}/build/doc/html";

    touch .nojekyll;

    git add --all .;

    git commit -m "Documentation build from Travis for commit ${TRAVIS_COMMIT}";

    git remote add upstream https://${GH_TOKEN}@github.com/${TRAVIS_REPO_SLUG};
    git push --quiet --force upstream gh-pages > /dev/null 2>&1;
}


cd "${TRAVIS_BUILD_DIR}";
gh_pages_prepare;
gh_pages_generate;
gh_pages_update;
