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

if [ "$TRAVIS_OS_NAME" = "linux" ]; then
	image=espressomd/espresso-$image:latest
	docker run $ci_env -u espresso --env-file $ENV_FILE -v ${PWD}:/travis -it $image /bin/bash -c "cp -r /travis .; cd travis && maintainer/travis/build_cmake.sh" || exit 1
elif [ "$TRAVIS_OS_NAME" = "osx" ]; then
	brew install cmake || brew upgrade cmake
	case "$image" in
		python3)
			brew install python3
			pip3 install h5py
			pip3 install cython
			pip3 install numpy
			pip3 install pep8
			PYTHON_VERSION=$(python3 --version 2>&1 | awk '{print $2}' | awk -F . '{print $1"."$2}')
			export cmake_params="-DPYTHON_EXECUTABLE=$(which python3) $cmake_params"
		;;
		*)
			brew install python
			pip install h5py
			pip install cython
			pip install numpy
			pip install pep8
		;;
	esac
	brew install boost-mpi
	brew install fftw
	travis_wait brew install hdf5 --with-mpi

	export TMPDIR=/tmp
	maintainer/travis/build_cmake.sh || exit 1
fi
