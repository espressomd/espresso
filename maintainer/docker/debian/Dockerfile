FROM debian:jessie
MAINTAINER Florian Weik <fweik@icp.uni-stuttgart.de>
ENV DEBIAN_FRONTEND noninteractive
RUN apt-get update
RUN apt-get install -y apt-utils
RUN apt-get install -y cmake 
RUN apt-get install -y build-essential
RUN apt-get install -y clang
RUN apt-get install -y openmpi-bin
RUN apt-get install -y libfftw3-dev
RUN apt-get install -y libboost-dev libboost-serialization-dev libboost-mpi-dev libboost-filesystem-dev libboost-test-dev
RUN apt-get install -y python python-numpy python-pyvtk
RUN apt-get install -y tcl-dev
RUN apt-get install -y git
RUN apt-get install -y pep8
RUN apt-get install -y python-pip
RUN apt-get install -y libpython-dev
RUN pip install cython

RUN useradd -m espresso
USER espresso
WORKDIR /home/espresso

