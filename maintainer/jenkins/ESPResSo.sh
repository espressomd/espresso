#!/bin/bash --login -e
source maintainer/jenkins/common.sh

if [ -n "$JENKINS_URL" ]; then
    echo "Pushing to local git clone..."
    git push --force /home/jenkins/espresso.git master
else
    echo "Would push to local git clone in Jenkins..."
fi


bootstrap

start "CONFIGURE"
./configure CPU_COUNT=4
end "CONFIGURE"

use_myconfig maxset

start "BUILD"
make -j 4
end "BUILD"

check
doc
dist

start "GITPUSH"
# push to savannah
if [ -z "$CHECK_UNSTABLE" ]; then
  echo "Pushing to Savannah..."
  if [ -n "$JENKINS_URL" ]; then
      git push ssh://olenz@git.sv.gnu.org/srv/git/espressomd.git master
  else
      echo "Would push to Savannah in Jenkins..."
  fi
else
  echo "Build is unstable, not pushing to Savannah."
fi
end "GITPUSH"
