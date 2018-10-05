#!/bin/bash

commits="c1d850ef3c4 c72cb35c3c5 eaa84cd1c92"
configs="myconfig-minimal.hpp ../src/core/myconfig-default.hpp ../maintainer/configs/maxset.hpp"

cd "$(git rev-parse --show-toplevel)"
mkdir -p build
cd build

cat > myconfig-minimal.hpp << EOF
#define ELECTROSTATICS
#define LENNARD_JONES
#define MASS
EOF

rm -f benchmarks.log
cat > benchmarks.csv << EOF
"commit","config","script","arguments","cores","MPI","mean","ci","steps_per_tick","duration","E1","E2","E3"
EOF

for commit in ${commits}
do
  rm -rf src/ maintainer/
  git checkout ${commit} ../src
  cmake .. -DWITH_BENCHMARKS=ON
  for config in ${configs}
  do
    cp ${config} myconfig.hpp
    make -j$(nproc)
    rm -f benchmarks.csv.part
    touch benchmarks.csv.part
    make benchmark 2>&1 | tee benchmarks.log
    sed -ri "s/^/\"${commit}\",\"$(basename ${config})\",/" benchmarks.csv.part
    cat benchmarks.csv.part >> benchmarks.csv
  done
done

git checkout HEAD ../src

