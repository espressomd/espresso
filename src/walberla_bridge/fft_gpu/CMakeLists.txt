#
# Copyright (C) 2020-2024 The ESPResSo project
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

add_library(espresso_fft_gpu OBJECT)
add_library(espresso::fft_gpu ALIAS espresso_fft_gpu)
target_link_libraries(
    espresso_fft_gpu
    PRIVATE
      espresso_cuda_flags
      Boost::mpi
      MPI::MPI_CXX
      Heffte::Heffte
      $<$<BOOL:${ESPRESSO_BUILD_WITH_SHARED_MEMORY_PARALLELISM}>:OpenMP::OpenMP_CXX>
      CUDA::cuda_driver
      CUDA::cudart
      walberla::core
      walberla::domain_decomposition
      walberla::blockforest
      walberla::communication
      walberla::field
      walberla::gpu
      espresso::walberla_codegen_cuda
      #espresso::walberla_cuda
      espresso_utils)
set_target_properties(espresso_fft_gpu PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries(espresso_walberla_cuda PRIVATE espresso_fft_gpu espresso_cuda_flags Heffte::Heffte)
target_sources(espresso_fft_gpu PRIVATE FFT_CUDA.cu)