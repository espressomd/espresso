#
# Copyright (C) 2020-2023 The ESPResSo project
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

include(unit_test)

function(ESPRESSO_WALBERLA_UNIT_TEST)
  cmake_parse_arguments(TEST "" "NAME;NUM_PROC" "SRC;DEPENDS" ${ARGN})
  unit_test(NAME ${TEST_NAME} NUM_PROC ${TEST_NUM_PROC} SRC ${TEST_SRC} DEPENDS
            ${TEST_DEPENDS} espresso::walberla espresso::walberla::cpp_flags
            espresso::utils)
  set_target_properties(${TEST_NAME} PROPERTIES CXX_CLANG_TIDY "")
  target_include_directories(${TEST_NAME} PRIVATE ${WALBERLA_INCLUDE_DIRS}
                                                  ${walberla_BINARY_DIR}/src)
  target_link_libraries(${TEST_NAME} PRIVATE ${WALBERLA_LIBS})
endfunction()

espresso_walberla_unit_test(NAME ResourceManager_test SRC
                            ResourceManager_test.cpp)

espresso_walberla_unit_test(NAME lb_kernels_unit_tests SRC
                            lb_kernels_unit_tests.cpp)

espresso_walberla_unit_test(NAME ek_kernels_unit_tests SRC
                            ek_kernels_unit_tests.cpp)

espresso_walberla_unit_test(
  NAME LatticeWalberla_unit_tests SRC LatticeWalberla_unit_tests.cpp DEPENDS
  Boost::boost Boost::mpi NUM_PROC 2)

espresso_walberla_unit_test(
  NAME LBWalberlaImpl_unit_tests SRC LBWalberlaImpl_unit_tests.cpp DEPENDS
  Boost::boost Boost::mpi NUM_PROC 2)

espresso_walberla_unit_test(
  NAME LBWalberlaImpl_bspline_tests SRC LBWalberlaImpl_bspline_tests.cpp
  DEPENDS Boost::mpi NUM_PROC 2)

if(NOT (ESPRESSO_BUILD_WITH_ASAN OR ESPRESSO_BUILD_WITH_UBSAN))
  espresso_walberla_unit_test(
    NAME LBWalberlaImpl_statistical_tests SRC
    LBWalberlaImpl_statistical_tests.cpp DEPENDS Boost::mpi)
endif()

espresso_walberla_unit_test(NAME LBWalberlaImpl_flow_tests SRC
                            LBWalberlaImpl_flow_tests.cpp DEPENDS Boost::mpi)

espresso_walberla_unit_test(NAME LBWalberlaImpl_lees_edwards_test SRC
                            LBWalberlaImpl_lees_edwards.cpp DEPENDS Boost::mpi)

espresso_walberla_unit_test(
  NAME EKinWalberlaImpl_unit_tests SRC EKinWalberlaImpl_unit_tests.cpp DEPENDS
  Boost::boost Boost::mpi NUM_PROC 2)
