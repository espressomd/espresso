whoami;
nvidia-smi;
cd ${CI_PROJECT_DIR};
cp maintainer/configs/maxset.hpp myconfig.hpp;
mkdir build && cd build;
cmake .. -DCUDA_NVCC_FLAGS="-D_FORCE_INLINES";
make -j 4 && make check_python;
result=$?;
if (( result == 0 )); then
    STATUS="success"
else
    STATUS="failure"
    cat testsuite/python/Testing/Temporary/LastTest.log
fi
GIT_COMMIT=$(git rev-parse HEAD)
URL=$(echo "https://gitlab.icp.uni-stuttgart.de/kai/espresso-build/pipelines/${CI_JOB_ID}"
curl "https://api.github.com/repos/kaiszuttor/espresso/statuses/$GIT_COMMIT?access_token=$GITHUB_TOKEN" \
      -H "Content-Type: application/json" \
      -X POST \
      -d "{\"state\": \"$STATUS\", \"description\": \"ICP build\", \"target_url\": \"$CI_ENVIRONMENT_URL\"}"
