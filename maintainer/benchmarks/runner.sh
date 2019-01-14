#!/bin/bash

cd "$(git rev-parse --show-toplevel)"
mkdir -p build
cd build

# manage headers files with different features
configs="myconfig-minimal.hpp myconfig-default.hpp myconfig-maxset.hpp"
cat > myconfig-minimal.hpp << EOF
#define ELECTROSTATICS
#define LENNARD_JONES
#define MASS
EOF
cp ../src/core/myconfig-default.hpp myconfig-default.hpp
sed 's/#define ADDITIONAL_CHECKS//' ../maintainer/configs/maxset.hpp > myconfig-maxset.hpp

# prepare build area
rm -rf src/ maintainer/
cmake -DWITH_BENCHMARKS=ON ..
cat > benchmarks.csv << EOF
"config","script","arguments","cores","MPI","mean","ci","steps_per_tick","duration","E1","E2","E3"
EOF

# run benchmarks
for config in ${configs}
do
  echo "### ${config}" >> benchmarks.log
  cp ${config} myconfig.hpp
  make -j$(nproc)
  rm -f benchmarks.csv.part
  touch benchmarks.csv.part
  make benchmark 2>&1 | tee -a benchmarks.log
  sed -ri "s/^/\"$(basename ${config})\",/" benchmarks.csv.part
  cat benchmarks.csv.part >> benchmarks.csv
done

rm benchmarks.csv.part

