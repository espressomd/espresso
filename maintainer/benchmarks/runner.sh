#!/bin/bash

cd "$(git rev-parse --show-toplevel)"
mkdir -p build
cd build

# manage configuration files with different features
configs="empty.hpp default.hpp maxset.hpp"
cp ../maintainer/configs/empty.hpp .
cp ../src/config/myconfig-default.hpp default.hpp
cp ../maintainer/configs/maxset.hpp .

# process configuration files
for config in ${configs}; do
  # add minimal features for the benchmarks to run
  echo -e "#define ELECTROSTATICS\n#define LENNARD_JONES\n#define MASS\n$(cat ${config})" > "${config}"
  # remove checks
  sed -ri "s/#define ADDITIONAL_CHECKS//" "${config}"
done

# prepare build area
rm -rf src/ maintainer/
cmake -DWITH_BENCHMARKS=ON -DTEST_TIMEOUT=600 -DWITH_CUDA=OFF -DWITH_CCACHE=OFF ..
cat > benchmarks.csv << EOF
"config","script","arguments","cores","MPI","mean","ci","steps_per_tick","duration"
EOF

# run benchmarks
for config in ${configs}; do
  echo "### ${config}" >> benchmarks.log
  cp ${config} myconfig.hpp
  make -j$(nproc)
  rm -f benchmarks.csv.part
  touch benchmarks.csv.part
  make -j$(nproc) benchmark 2>&1 | tee -a benchmarks.log
  sed -ri "s/^/\"$(basename ${config})\",/" benchmarks.csv.part
  cat benchmarks.csv.part >> benchmarks.csv
done

rm benchmarks.csv.part

