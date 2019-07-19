#!/bin/bash

abort() {
  echo "An error occurred in runner.sh, exiting now" >&2
  echo "Command that failed: ${BASH_COMMAND}" >&2
  exit 1
}

trap abort EXIT
set -e

# manage configuration files with different features
configs="empty.hpp default.hpp maxset.hpp"
cp ../maintainer/configs/empty.hpp .
cp ../src/config/myconfig-default.hpp default.hpp
cp ../maintainer/configs/maxset.hpp .

# process configuration files
for config in ${configs}; do
  # add minimal features for the benchmarks to run
  sed -i '1 i\#define ELECTROSTATICS\n#define LENNARD_JONES\n#define MASS\n' "${config}"
  # remove checks
  sed -ri "s/#define ADDITIONAL_CHECKS//" "${config}"
done

cat > benchmarks.csv << EOF
"config","script","arguments","cores","MPI","mean","ci","nsteps","duration"
EOF

# run benchmarks
for config in ${configs}; do
  echo "### ${config}" | tee -a benchmarks.log
  cp ${config} myconfig.hpp
  rm -rf src/ maintainer/
  cmake -DWITH_BENCHMARKS=ON -DTEST_TIMEOUT=600 -DWITH_CUDA=OFF -DWITH_CCACHE=OFF ..
  make -j$(nproc)
  rm -f benchmarks.csv.part
  touch benchmarks.csv.part
  make benchmark 2>&1 | tee -a benchmarks.log
  sed -ri "s/^/\"$(basename ${config})\",/" benchmarks.csv.part
  cat benchmarks.csv.part >> benchmarks.csv
done

rm benchmarks.csv.part

trap : EXIT
