#!/bin/bash

ENV_FILE=$(mktemp esXXXXXXX.env)

cat > $ENV_FILE <<EOF
insource=$insource
configure_params=$configure_params
with_fftw=$with_fftw
with_tcl=$with_tcl
with_python_interface=$with_python_interface
myconfig=$myconfig
check_procs=$check_procs
make_check=$make_check
EOF

if [ -z "$image" ]; then
	image=ubuntu
fi

docker run -u espresso --env-file $ENV_FILE -v ${PWD}:/travis -it espressomd/buildenv-espresso-${image} /bin/bash -c "git clone /travis && cd travis && maintainer/travis/build_cmake.sh"
