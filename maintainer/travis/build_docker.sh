#!/bin/bash

ENV_FILE=$(mktemp esXXXXXXX.env)

cat > $ENV_FILE <<EOF
insource=$insource
cmake_params=$cmake_params
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

if [ "$TRAVIS_OS_NAME" = "linux" ]; then
	image=espressomd/buildenv-espresso-$image
	docker run -u espresso --env-file $ENV_FILE -v ${PWD}:/travis -it $image /bin/bash -c "git clone /travis && cd travis && maintainer/travis/build_cmake.sh"
elif [ "$TRAVIS_OS_NAME" = "osx" ]; then
	brew install cmake || brew upgrade cmake
	case "$image" in
		python3)
			brew install python3
			pip3 install cython
			pip3 install numpy
			pip3 install pep8
			PYTHON_VERSION=$(python3 --version 2>&1 | awk '{print $2}' | awk -F . '{print $1"."$2}')
			export cmake_params="-DPYTHON_EXECUTABLE=$(which python3) $cmake_params"
		;;
		*)
			brew install python
			pip install cython
			pip install numpy
			pip install pep8
		;;
	esac
	brew install boost-mpi
	brew install fftw

	maintainer/travis/build_cmake.sh
fi
