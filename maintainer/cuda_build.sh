#!/usr/bin/env bash

GIT_COMMIT=$(git rev-parse HEAD)
URL=$(echo "https://gitlab.icp.uni-stuttgart.de/espressomd/espresso/-/jobs/${CI_JOB_ID}")
STATUS="pending"
curl "https://api.github.com/repos/espressomd/espresso/statuses/$GIT_COMMIT?access_token=$GITHUB_TOKEN" \
      -H "Content-Type: application/json" \
      -X POST \
      -d "{\"state\": \"$STATUS\", \"context\": \"ICP CUDA build\", \"target_url\": \"$URL\"}"

nvidia-smi
cd ${CI_PROJECT_DIR}
cp maintainer/configs/maxset.hpp myconfig.hpp
mkdir build && cd build
cmake .. -DCUDA_NVCC_FLAGS="-D_FORCE_INLINES"
make -j 4 && make check_python
result=$?
echo $result
if (( $result == 0 )); then
    STATUS="success"
else
    STATUS="failure"
    cat testsuite/python/Testing/Temporary/LastTest.log
fi
echo "result was $result\n"
GIT_COMMIT=$(git rev-parse HEAD)
URL=$(echo "https://gitlab.icp.uni-stuttgart.de/espressomd/espresso/-/jobs/${CI_JOB_ID}")
curl "https://api.github.com/repos/espressomd/espresso/statuses/$GIT_COMMIT?access_token=$GITHUB_TOKEN" \
      -H "Content-Type: application/json" \
      -X POST \
      -d "{\"state\": \"$STATUS\", \"context\": \"ICP CUDA build\", \"target_url\": \"$URL\"}"
exit $result
