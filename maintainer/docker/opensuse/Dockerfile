FROM opensuse:42.3
MAINTAINER Florian Weik <fweik@icp.uni-stuttgart.de>
RUN zypper -n --gpg-auto-import-keys refresh \
&& zypper -n --gpg-auto-import-keys install \
  patterns-openSUSE-devel_C_C++ gcc-c++ \
  cmake \
  openmpi-devel Modules \
  fftw3-devel \
  boost-devel libboost_serialization1_54_0 libboost_mpi1_54_0  libboost_filesystem1_54_0 libboost_test1_54_0 \
  python python-numpy python-numpy-devel python-h5py \
  git \
  python-pep8 \
  python-pip \
  python-Cython \
  python-devel \
  hdf5-openmpi-devel libhdf5-10-openmpi
RUN useradd -m espresso
USER espresso
WORKDIR /home/espresso
