#
# Copyright (C) 2022-2026 The ESPResSo project
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

add_subdirectory(generated_kernels)

target_sources(espresso_walberla PRIVATE lb_walberla_init.cpp)
set_property(
  SOURCE lb_walberla_init.cpp TARGET_DIRECTORY espresso_walberla APPEND
  PROPERTY COMPILE_OPTIONS $<$<COMPILE_LANG_AND_ID:CXX,GNU>:-fno-analyzer>)

if(ESPRESSO_BUILD_WITH_CUDA)
  target_sources(espresso_walberla PRIVATE lb_walberla_init.cu)
  set_property(
    SOURCE lb_walberla_init.cu TARGET_DIRECTORY espresso_walberla
    APPEND
    PROPERTY
      COMPILE_OPTIONS
      $<$<AND:$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>,$<STREQUAL:${CMAKE_CUDA_HOST_COMPILER_ID},GNU>>:-Xcompiler=-fno-analyzer>
  )
endif()
