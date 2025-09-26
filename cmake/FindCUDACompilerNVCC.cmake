#
# Copyright (C) 2009-2025 The ESPResSo project
# Copyright (C) 2009,2010
#   Max-Planck-Institute for Polymer Research, Theory Group
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
#

# include the toolkit libraries and declare a custom
# `add_library()` wrapper function named `espresso_add_gpu_library()`.

set(CUDA_LINK_LIBRARIES_KEYWORD PUBLIC)
set(CUDA_PROPAGATE_HOST_FLAGS OFF)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  CUDACompilerNVCC REQUIRED_VARS CMAKE_CUDA_COMPILER VERSION_VAR
  CMAKE_CUDA_COMPILER_VERSION)
