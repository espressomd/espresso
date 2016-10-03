FROM nvidia/cuda:7.5-devel-ubuntu14.04
MAINTAINER Florian Weik <fweik@icp.uni-stuttgart.de>
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update && apt-get install -y \
	apt-utils \
	build-essential \
	curl \
	openmpi-bin \
	libfftw3-dev \
	libboost-dev libboost-serialization-dev libboost-mpi-dev libboost-filesystem-dev libboost-test-dev \
	python python-numpy \
	tcl-dev \
	git \
	pep8 \
	python-pyvtk \
	python-pip \
	libpython-dev \
&& pip install cython \
&& apt-get clean \
&& rm -rf /var/lib/apt/lists/*

RUN cd /usr/local \
&& curl -sL https://cmake.org/files/v3.1/cmake-3.1.0-Linux-x86_64.tar.gz | tar --strip-components=1 -xz

RUN useradd -m espresso
USER espresso
WORKDIR /home/espresso

