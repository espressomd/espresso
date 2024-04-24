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

add_library(espresso_walberla_cpp_flags INTERFACE)
add_library(espresso::walberla::cpp_flags ALIAS espresso_walberla_cpp_flags)
target_link_libraries(
  espresso_walberla_cpp_flags
  INTERFACE espresso::cpp_flags
            $<$<BOOL:${ESPRESSO_BUILD_WITH_WALBERLA_AVX}>:espresso::avx_flags>)

function(espresso_configure_walberla_target)
  set(TARGET_NAME ${ARGV0})
  set_target_properties(${TARGET_NAME} PROPERTIES CXX_CLANG_TIDY "")
  target_link_libraries(${TARGET_NAME} PRIVATE ${WALBERLA_LIBS})
  target_include_directories(
    ${TARGET_NAME} PUBLIC include PRIVATE ${WALBERLA_INCLUDE_DIRS}
                                          ${walberla_BINARY_DIR}/src)
  install(TARGETS ${TARGET_NAME}
          LIBRARY DESTINATION ${ESPRESSO_INSTALL_PYTHON}/espressomd)
endfunction()

add_library(espresso_walberla SHARED)
add_library(espresso::walberla ALIAS espresso_walberla)

espresso_configure_walberla_target(espresso_walberla)

target_link_libraries(espresso_walberla PUBLIC MPI::MPI_CXX espresso::utils
                      PRIVATE espresso::walberla::cpp_flags)

if(ESPRESSO_BUILD_WITH_CUDA AND WALBERLA_BUILD_WITH_CUDA)
  espresso_add_gpu_library(espresso_walberla_cuda SHARED)
  add_library(espresso::walberla_cuda ALIAS espresso_walberla_cuda)
  espresso_configure_walberla_target(espresso_walberla_cuda)
  target_link_libraries(
    espresso_walberla_cuda PUBLIC espresso::utils espresso::walberla
    PRIVATE CUDA::cuda_driver CUDA::cudart)
endif()

add_subdirectory(src)
if(ESPRESSO_BUILD_TESTS)
  add_subdirectory(tests)
endif()
