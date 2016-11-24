#!/bin/bash

ENV_FILE=$(mktemp esXXXXXXX.env)

cat > $ENV_FILE <<EOF
insource=$insource
cmake_params=$cmake_params
with_fftw=$with_fftw
with_python_interface=yes
myconfig=$myconfig
check_procs=$check_procs
make_check=$make_check
EOF

if [ -z "$image" ]; then
	image=ubuntu
fi

if [ "$TRAVIS_OS_NAME" = "linux" ]; then
	image=espressomd/buildenv-espresso-$image:python
	docker run -u espresso --env-file $ENV_FILE -v ${PWD}:/travis -it $image /bin/bash -c "git clone --recursive /travis && cd travis && maintainer/travis/build_cmake.sh"
elif [ "$TRAVIS_OS_NAME" = "osx" ]; then
#	brew update
#	brew upgrade
	brew install cmake
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
	brew install openmpi
	brew install fftw
    #brew install homebrew/science/hdf5 --with-mpi
    HDF5_VERSION=1.8.17
    export HDF5_ROOT=$HOME/hdf5
    wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-${HDF5_VERSION}/src/hdf5-${HDF5_VERSION}.tar.gz && tar xfvz hdf5-${HDF5_VERSION}.tar.gz &&
    cd hdf5-${HDF5_VERSION} && mkdir build && cd build && cmake .. -DHDF5_ENABLE_PARALLEL=ON -DHDF5_BUILD_CPP_LIB=OFF -DCMAKE_INSTALL_PREFIX=$HOME/hdf5 &&
    make -j4 &> /dev/null && make install && cd ../../ ;
    # The binary version of Boost comes without MPI support, so we have to
    # compile it ourselves.
	# Boost takes a long time to install, so we have Travis-CI cache it.
	BOOST_VERSION=$(brew info --json=v1 boost | python -m json.tool | grep linked_keg | awk -F '"' '{print $4}')
	brew unlink boost
	brew info --json=v1 boost | python -m json.tool | grep -B 4 "\"version\": \"$BOOST_VERSION\"" | grep -q 'with-mpi' || brew uninstall --ignore-dependencies boost
	brew install boost --without-single --with-mpi &
	BOOST_PID=$!
	START_TIME=$(date +%s)
	while kill -0 $BOOST_PID; do
		# boost takes a while to compile, during which there 
		# is no output, causing Travis-CI to abort eventually
		echo "Boost is being built..."
		sleep 30
	done
	brew link boost
	if [ "$(echo $(date +%s)-$START_TIME | bc)" -gt 600 ]; then
		# Boost was built, i.e. it took more than 10 minutes, so
		# abort build because we won't be able to complete it
		# within the time limit
		echo
		echo "============================================================================="
		echo "Boost was built from source, but this took so long that we won't have time to"
		echo "finish the actual build process before the Travis timeout kills it.
		echo "Please manually retry the current Travis build!
		echo "============================================================================="
		echo
		exit 2
	fi

	maintainer/travis/build_cmake.sh
fi
