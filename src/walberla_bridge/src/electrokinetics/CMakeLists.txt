#
# Copyright (C) 2022-2023 The ESPResSo project
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
add_subdirectory(reactions)

target_sources(espresso_walberla PRIVATE ek_walberla_init.cpp)
target_sources(espresso_walberla PRIVATE ek_poisson_none_init.cpp)
if(ESPRESSO_BUILD_WITH_WALBERLA_FFT)
  target_sources(espresso_walberla PRIVATE ek_poisson_fft_init.cpp)
endif()
if(ESPRESSO_BUILD_WITH_CUDA AND WALBERLA_BUILD_WITH_CUDA)
  if(ESPRESSO_BUILD_WITH_WALBERLA_FFT)
    add_library(espresso_fft_gpu OBJECT)
    add_library(espresso::fft_gpu ALIAS espresso_fft_gpu)
    target_sources(espresso_fft_gpu PRIVATE FFT_CUDA.cu)
    target_link_libraries(
      espresso_fft_gpu PRIVATE espresso_cuda_flags Heffte::Heffte walberla::gpu
                               espresso::walberla_codegen_cuda espresso_utils)
    set_target_properties(espresso_fft_gpu PROPERTIES POSITION_INDEPENDENT_CODE
                                                      ON)
    target_link_libraries(espresso_walberla_cuda PRIVATE espresso_fft_gpu)
    target_sources(espresso_walberla_cuda PRIVATE ek_poisson_fft_gpu_init.cu)
  endif()
  target_sources(espresso_walberla_cuda PRIVATE ek_walberla_init.cu)
endif()
