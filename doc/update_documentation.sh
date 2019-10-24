#!/bin/sh

set -e;

# Gather coverage data (always has to be done for CodeCov)
lcov --directory . --capture --output-file coverage.info;
lcov --remove coverage.info '/usr/*' '*/examples/*' --output-file coverage.info;

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


# Build the docs
make -k -j 2 doxygen;

# Generate HTML coverage report
genhtml coverage.info --output-directory doc/coverage/html/;

cd "doc/";

touch .nojekyll;

git add --all ".nojekyll" "doxygen/html" "coverage/html";

git commit --quiet -m "Documentation build from Travis for commit ${TRAVIS_COMMIT}";

git remote add upstream https://${GH_TOKEN}@github.com/${TRAVIS_REPO_SLUG};
git push --quiet --force upstream gh-pages > /dev/null 2>&1;
