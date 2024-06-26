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

include(espresso_unit_test)

function(ESPRESSO_ADD_TEST)
  cmake_parse_arguments(TEST "" "SRC;NAME;NUM_PROC" "DEPENDS" ${ARGN})
  espresso_unit_test(
    SRC ${TEST_SRC} NAME ${TEST_NAME} NUM_PROC ${TEST_NUM_PROC} DEPENDS
    ${TEST_DEPENDS} espresso::walberla espresso::utils)
  if(${TEST_SRC} MATCHES ".*\.cu$")
    target_link_libraries(${TEST_NAME} PRIVATE espresso::walberla::cuda_flags
                                               espresso::walberla_cuda)
  else()
    target_link_libraries(${TEST_NAME} PRIVATE espresso::walberla::cpp_flags)
  endif()
  set_target_properties(${TEST_NAME} PROPERTIES CXX_CLANG_TIDY "")
  target_include_directories(${TEST_NAME} PRIVATE ${WALBERLA_INCLUDE_DIRS}
                                                  ${walberla_BINARY_DIR}/src)
  target_link_libraries(${TEST_NAME} PRIVATE ${WALBERLA_LIBS})
endfunction()

espresso_add_test(SRC ResourceManager_test.cpp)
espresso_add_test(SRC lb_kernels_unit_tests.cpp)
espresso_add_test(SRC ek_kernels_unit_tests.cpp)
espresso_add_test(SRC LatticeWalberla_unit_tests.cpp DEPENDS Boost::mpi
                  NUM_PROC 2)
espresso_add_test(SRC LBWalberlaImpl_unit_tests.cpp DEPENDS Boost::mpi NUM_PROC
                  2)
espresso_add_test(SRC LBWalberlaImpl_bspline_tests.cpp DEPENDS Boost::mpi
                  NUM_PROC 2)
espresso_add_test(SRC LBWalberlaImpl_flow_tests.cpp DEPENDS Boost::mpi)
espresso_add_test(SRC LBWalberlaImpl_lees_edwards_tests.cpp DEPENDS Boost::mpi)
espresso_add_test(SRC EKinWalberlaImpl_unit_tests.cpp DEPENDS Boost::mpi
                  NUM_PROC 2)

if(NOT (ESPRESSO_BUILD_WITH_ASAN OR ESPRESSO_BUILD_WITH_UBSAN))
  espresso_add_test(SRC LBWalberlaImpl_statistical_tests.cpp DEPENDS Boost::mpi)
endif()

if(ESPRESSO_BUILD_WITH_CUDA AND WALBERLA_BUILD_WITH_CUDA)
  espresso_add_test(SRC LBWalberlaImpl_field_accessors_tests.cu DEPENDS
                    Boost::mpi)
endif()
