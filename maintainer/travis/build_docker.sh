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
EOF

docker run -u espresso --env-file $ENV_FILE -v ${PWD}:/travis -it espressomd/travis-build /bin/bash -c "git clone /travis && cd travis && maintainer/travis/build_cmake.sh"
