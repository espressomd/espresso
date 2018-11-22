#!/bin/bash

# list of commits to benchmark
commits="c1d850ef3c4 c72cb35c3c5 eaa84cd1c92"

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
# disable features interfering with LJ benchmarks
sed -ri 's/#define +(THOLE|COLLISION_DETECTION|GHOSTS_HAVE_BONDS)/\/\/#define \1/' myconfig-*.hpp

# output files
rm -f benchmarks.log
cat > benchmarks.csv << EOF
"commit","config","script","arguments","cores","MPI","mean","ci","steps_per_tick","duration","E1","E2","E3"
EOF

# run benchmarks
for commit in ${commits}
do
  rm -rf src/ maintainer/
  git checkout ${commit} ../src ../libs
  cmake .. -DWITH_BENCHMARKS=ON
  for config in ${configs}
  do
    echo "### ${commit} ${config}" >> benchmarks.log
    cp ${config} myconfig.hpp
    make -j$(nproc)
    rm -f benchmarks.csv.part
    touch benchmarks.csv.part
    make benchmark 2>&1 | tee -a benchmarks.log
    sed -ri "s/^/\"${commit}\",\"$(basename ${config})\",/" benchmarks.csv.part
    cat benchmarks.csv.part >> benchmarks.csv
  done
done

# restore files
git checkout HEAD ../src ../libs

