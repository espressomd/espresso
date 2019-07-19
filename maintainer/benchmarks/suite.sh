#!/bin/bash

# list of commits to benchmark
commits="HEAD"

# list of directories to checkout
directories="../src ../libs"

abort() {
    echo "An error occurred in suite.sh, exiting now" >&2
    echo "Command that failed: ${BASH_COMMAND}" >&2
    cleanup
    exit 1
}

cleanup() {
    # restore files
    git checkout HEAD ${directories}
}

trap abort EXIT
set -e

# move to top-level directory
cd "$(git rev-parse --show-toplevel)"
build_dir="$(realpath build-benchmarks)"

# move to build directory
if [ -d "${build_dir}" ]; then
  rm -rf "${build_dir}"
fi
mkdir "${build_dir}"
cd "${build_dir}"

# prepare output files
rm -f benchmarks.log
cat > benchmarks_suite.csv << EOF
"commit","config","script","arguments","cores","MPI","mean","ci","steps_per_tick","duration"
EOF

# run benchmarks
for commit in ${commits}
do
  echo "### commit ${commit}" >> benchmarks.log
  git checkout ${commit} ${directories}
  bash ../maintainer/benchmarks/runner.sh
  sed -ri "s/^/\"${commit}\",/" benchmarks.csv
  tail -n +2 benchmarks.csv >> benchmarks_suite.csv
done

rm benchmarks.csv

trap : EXIT
cleanup
