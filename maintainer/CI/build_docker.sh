#!/usr/bin/env sh
# Copyright (C) 2016-2019 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

ENV_FILE=$(mktemp esXXXXXXX.env)

cat > "${ENV_FILE}" <<EOF
cmake_params=${cmake_params}
with_fftw=${with_fftw}
with_coverage=false
with_cuda=true
CC=gcc-8
CXX=g++-8
check_procs=${check_procs}
make_check=${make_check}
build_type=RelWithAssert
myconfig=maxset
EOF
image="espressomd/docker-ubuntu-20.04:06b6216c7aa3555bcf28c90734dbb84e7285c96f"

docker run -u espresso --env-file "${ENV_FILE}" -v "${PWD}:/travis" -it "${image}" /bin/bash -c "cp -r /travis .; cd travis && maintainer/CI/build_cmake.sh" || exit 1
