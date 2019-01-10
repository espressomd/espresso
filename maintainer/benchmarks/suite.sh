#!/bin/bash

# list of commits to benchmark
commits="HEAD"

cd "$(git rev-parse --show-toplevel)"
mkdir -p build
cd build

# prepare output files
rm -f benchmarks.log
cat > benchmarks_suite.csv << EOF
"commit","config","script","arguments","cores","MPI","mean","ci","steps_per_tick","duration","E1","E2","E3"
EOF

# run benchmarks
for commit in ${commits}
do
  echo "### commit ${commit}" >> benchmarks.log
  git checkout ${commit} ../src ../libs
  bash ../maintainer/benchmarks/runner.sh
  sed -ri "s/^/\"${commit}\",/" benchmarks.csv
  tail -n +2 benchmarks.csv >> benchmarks_suite.csv
done

rm benchmarks.csv

# restore files
git checkout HEAD ../src ../libs

