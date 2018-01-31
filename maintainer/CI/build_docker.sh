#!/bin/bash

ENV_FILE=$(mktemp esXXXXXXX.env)

ci_env=""
if [[ ! -z ${with_coverage+x} ]]; then 
  ci_env=`bash <(curl -s https://codecov.io/env)`
fi

cat > $ENV_FILE <<EOF
insource=$insource
cmake_params=$cmake_params
with_fftw=$with_fftw
with_python_interface=yes
with_coverage=$with_coverage
myconfig=$myconfig
check_procs=$check_procs
make_check=$make_check
EOF

if [ -z "$image" ]; then
	image=ubuntu
fi

image=espressomd/espresso-$image:latest
docker run $ci_env -u espresso --env-file $ENV_FILE -v ${PWD}:/travis -it $image /bin/bash -c "cp -r /travis .; cd travis && maintainer/CI/build_cmake.sh" || exit 1
