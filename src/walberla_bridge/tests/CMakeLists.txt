#
# Copyright (C) 2020-2026 The ESPResSo project
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
  cmake_parse_arguments(TEST "" "SRC;NAME;NUM_PROC;NUM_THREADS" "DEPENDS"
                        ${ARGN})
  espresso_unit_test_make_name(${ARGV})
  espresso_unit_test(
    SRC
    ${TEST_SRC}
    NAME
    ${TEST_NAME}
    NUM_PROC
    ${TEST_NUM_PROC}
    NUM_THREADS
    ${TEST_NUM_THREADS}
    DEPENDS
    ${TEST_DEPENDS}
    $<$<BOOL:${ESPRESSO_BUILD_WITH_WALBERLA_AVX}>:espresso::avx_flags>
    $<$<BOOL:${ESPRESSO_BUILD_WITH_SHARED_MEMORY_PARALLELISM}>:OpenMP::OpenMP_CXX>
    espresso::walberla
    espresso::walberla_deps
    espresso::compiler_flags
    espresso::walberla::compiler_flags
    espresso::utils)
  target_include_directories(${TEST_NAME} PRIVATE ${WALBERLA_INCLUDE_DIRS}
                                                  ${walberla_BINARY_DIR}/src)
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
espresso_add_test(
  SRC PoissonSolver_test.cpp DEPENDS Boost::mpi
  $<$<BOOL:${ESPRESSO_BUILD_WITH_FFTW}>:Heffte::Heffte> NUM_PROC 2)

if(NOT (ESPRESSO_BUILD_WITH_ASAN OR ESPRESSO_BUILD_WITH_UBSAN))
  espresso_add_test(SRC LBWalberlaImpl_statistical_tests.cpp DEPENDS Boost::mpi)
endif()

espresso_add_test(SRC LBWalberlaImpl_field_accessors_tests.cpp DEPENDS
                  Boost::mpi)
if(ESPRESSO_BUILD_WITH_CUDA AND WALBERLA_BUILD_WITH_CUDA)
  set_source_files_properties(LBWalberlaImpl_field_accessors_tests.cpp
                              PROPERTIES LANGUAGE CUDA)
  espresso_add_test(SRC LBWalberlaImpl_field_accessors_tests.cpp SUFFIX cuda
                    DEPENDS Boost::mpi CUDA::cudart)
endif()
