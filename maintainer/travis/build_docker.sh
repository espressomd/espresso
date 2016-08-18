#!/bin/bash

docker run -u espresso -e insource=$insource -e configure_params=$configure_params -e with_fftw=$with_fftw -e with_tcl=$with_tcl -e with_python_interface=$with_python_interface -e myconfig=$myconfig -e check_procs=$check_procs -v ${PWD}:/travis -it espressomd/travis-build /bin/bash -c "git clone /travis && cd travis && maintainer/travis/build_cmake.sh"


