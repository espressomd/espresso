FROM opensuse:42.1
MAINTAINER Florian Weik <fweik@icp.uni-stuttgart.de>
RUN zypper -n --gpg-auto-import-keys refresh \
&& zypper -n --gpg-auto-import-keys install \
  patterns-openSUSE-devel_C_C++ \
  cmake \
  openmpi-devel \
  fftw3-devel \
  boost-devel libboost_serialization1_54_0 libboost_mpi1_54_0  libboost_filesystem1_54_0 libboost_test1_54_0 \
  python python-numpy python-PyVTK \
  tcl-devel \
  git \
  python-pep8 \
  python-pip \
  libpython-devel

RUN useradd -m espresso
USER espresso
WORKDIR /home/espresso
