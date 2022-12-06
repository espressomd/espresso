#
# Copyright (C) 2020-2022 The ESPResSo project
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

add_library(espresso_walberla SHARED)
add_library(espresso::walberla ALIAS espresso_walberla)

set_target_properties(espresso_walberla PROPERTIES CXX_CLANG_TIDY "")
target_include_directories(espresso_walberla PUBLIC include)

add_library(espresso_walberla_cpp_flags INTERFACE)
set_target_properties(espresso_walberla_cpp_flags PROPERTIES CXX_CLANG_TIDY "")
add_library(espresso::walberla::cpp_flags ALIAS espresso_walberla_cpp_flags)
if(ESPRESSO_BUILD_WITH_WALBERLA_USE_AVX)
  target_link_libraries(espresso_walberla_cpp_flags
                        INTERFACE espresso::avx_flags)
endif()

#if(ESPRESSO_BUILD_WITH_CUDA AND WALBERLA_BUILD_WITH_CUDA)
#  espresso_add_gpu_library(espresso_walberla_cuda SHARED
#                  ${WalberlaBridgeCodegenCuda_SRC})
#  add_library(espresso::walberla_cuda ALIAS espresso_walberla_cuda)
#  target_link_libraries(espresso_walberla PUBLIC espresso::walberla_cuda)
#endif()

target_link_libraries(
  espresso_walberla PUBLIC MPI::MPI_CXX espresso::utils
  PRIVATE Boost::boost espresso::cpp_flags espresso::walberla::cpp_flags
          ${WALBERLA_LIBS})
# TODO WALBERLA: break dependency on walberla headers in public API
target_include_directories(
  espresso_walberla PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${WALBERLA_INCLUDE_DIRS}
  PRIVATE ${walberla_BINARY_DIR}/src)

add_subdirectory(src)
# TODO WALBERLA: find why `TARGET_DIRECTORY espresso_walberla` doesn't work
# in subfolders where these source files are added to `espresso_walberla`
set_source_files_properties(
  src/electrokinetics/ek_walberla_init.cpp
  PROPERTIES COMPILE_FLAGS -Wno-unused-variable)

install(TARGETS espresso_walberla
        LIBRARY DESTINATION ${ESPRESSO_INSTALL_PYTHON}/espressomd)

if(ESPRESSO_BUILD_TESTS)
  set(CMAKE_CXX_CLANG_TIDY "")
  add_subdirectory(tests)
  set(CMAKE_CXX_CLANG_TIDY "${ESPRESSO_DEFAULT_CXX_CLANG_TIDY}")
endif()
