#!/usr/bin/env bash

GIT_COMMIT=$(git rev-parse HEAD)
URL=$(echo "https://gitlab.icp.uni-stuttgart.de/espressomd/espresso/-/jobs/${CI_JOB_ID}")
STATUS="pending"
curl "https://api.github.com/repos/espressomd/espresso/statuses/$GIT_COMMIT?access_token=$GITHUB_TOKEN" \
      -H "Content-Type: application/json" \
      -X POST \
      -d "{\"state\": \"$STATUS\", \"context\": \"ICP CUDA build\", \"target_url\": \"$URL\"}"

nvidia-smi
srcdir="${CI_PROJECT_DIR}" cmake_params='-DCUDA_NVCC_FLAGS="-D_FORCE_INLINES"' with_coverage="true" build_procs=4 maintainer/travis/build_cmake.sh
result=$?
echo $result
if (( $result == 0 )); then
    STATUS="success"
else
    STATUS="failure"
fi
echo "result was $result\n"
GIT_COMMIT=$(git rev-parse HEAD)
URL=$(echo "https://gitlab.icp.uni-stuttgart.de/espressomd/espresso/-/jobs/${CI_JOB_ID}")
curl "https://api.github.com/repos/espressomd/espresso/statuses/$GIT_COMMIT?access_token=$GITHUB_TOKEN" \
      -H "Content-Type: application/json" \
      -X POST \
      -d "{\"state\": \"$STATUS\", \"context\": \"ICP CUDA build\", \"target_url\": \"$URL\"}"
exit $result
