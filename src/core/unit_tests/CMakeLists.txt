#
# Copyright (C) 2010-2026 The ESPResSo project
# Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
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

include(espresso_unit_test)

espresso_unit_test(SRC RuntimeError_test.cpp DEPENDS Boost::serialization)
espresso_unit_test(SRC RuntimeErrorCollector_test.cpp DEPENDS espresso::core
                   Kokkos::kokkos Boost::mpi MPI::MPI_CXX NUM_PROC 2)
espresso_unit_test(SRC RuntimeErrorStream_test.cpp DEPENDS espresso::core
                   Kokkos::kokkos Boost::mpi MPI::MPI_CXX NUM_PROC 2)
# create a test matrix for the stand-alone unit test
espresso_unit_test_executable(
  NAME EspressoSystemStandAlone_test_matrix SRC
  EspressoSystemStandAlone_test.cpp DEPENDS espresso::core Kokkos::kokkos
  Boost::mpi MPI::MPI_CXX)
foreach(TEST_NUM_PROC 1 3 4)
  foreach(TEST_NUM_THREADS 1 2)
    if(${TEST_NUM_THREADS} LESS_EQUAL ${ESPRESSO_TEST_NT})
      if(${TEST_NUM_PROC} LESS_EQUAL ${ESPRESSO_TEST_NP})
        set(TEST_NAME
            EspressoSystemStandAlone_test_${TEST_NUM_PROC}_mpi_ranks_${TEST_NUM_THREADS}_omp_threads
        )
        espresso_unit_test_register(
          NAME ${TEST_NAME} TARGET EspressoSystemStandAlone_test_matrix
          NUM_PROC ${TEST_NUM_PROC} NUM_THREADS ${TEST_NUM_THREADS})
      endif()
    endif()
  endforeach()
endforeach()
espresso_unit_test(SRC EspressoSystem_test.cpp DEPENDS espresso::core
                   Kokkos::kokkos Boost::mpi NUM_PROC 2)
espresso_unit_test(SRC ResourceCleanup_test.cpp DEPENDS espresso::core
                   Kokkos::kokkos Boost::mpi NUM_PROC 2)
espresso_unit_test(SRC MpiCallbacks_test.cpp DEPENDS espresso::utils Boost::mpi
                   MPI::MPI_CXX NUM_PROC 2)
espresso_unit_test(SRC ParticleIterator_test.cpp DEPENDS espresso::utils)
espresso_unit_test(SRC link_cell_test.cpp DEPENDS espresso::utils)
espresso_unit_test(SRC Particle_test.cpp DEPENDS espresso::utils
                   Boost::serialization)
espresso_unit_test(SRC Particle_serialization_test.cpp DEPENDS espresso::utils
                   Boost::serialization)
espresso_unit_test(SRC rotation_test.cpp DEPENDS espresso::utils espresso::core
                   Kokkos::kokkos)
espresso_unit_test(SRC field_coupling_couplings_test.cpp DEPENDS
                   espresso::utils)
espresso_unit_test(SRC field_coupling_fields_test.cpp DEPENDS espresso::utils)
espresso_unit_test(SRC field_coupling_force_field_test.cpp DEPENDS
                   espresso::utils)
espresso_unit_test(SRC field_layout_helpers_test.cpp DEPENDS espresso::utils
                   espresso::core Kokkos::kokkos)
espresso_unit_test(SRC particle_reduction_test.cpp DEPENDS espresso::utils
                   espresso::core Kokkos::kokkos Boost::mpi)
espresso_unit_test(SRC dipolar_direct_sum_kernels_test.cpp DEPENDS
                   espresso::core espresso::utils Kokkos::kokkos)
espresso_unit_test(SRC periodic_fold_test.cpp)
espresso_unit_test(SRC grid_test.cpp DEPENDS espresso::core Kokkos::kokkos)
espresso_unit_test(SRC lees_edwards_test.cpp DEPENDS espresso::core
                   Kokkos::kokkos)
espresso_unit_test(SRC BoxGeometry_test.cpp DEPENDS espresso::core
                   Kokkos::kokkos)
espresso_unit_test(SRC LocalBox_test.cpp DEPENDS espresso::core Kokkos::kokkos)
espresso_unit_test(SRC Verlet_list_test.cpp DEPENDS espresso::core
                   Kokkos::kokkos NUM_PROC 4)
espresso_unit_test(SRC VerletCriterion_test.cpp DEPENDS espresso::core
                   Kokkos::kokkos)
espresso_unit_test(SRC thermostats_test.cpp DEPENDS espresso::core
                   Kokkos::kokkos)
espresso_unit_test(SRC random_test.cpp DEPENDS espresso::utils Kokkos::kokkos
                   Random123)
espresso_unit_test(SRC BondList_test.cpp DEPENDS espresso::core Kokkos::kokkos)
espresso_unit_test(SRC energy_test.cpp DEPENDS espresso::core Kokkos::kokkos)
espresso_unit_test(SRC bonded_interactions_map_test.cpp DEPENDS espresso::core
                   Kokkos::kokkos)
espresso_unit_test(SRC bond_breakage_test.cpp DEPENDS espresso::core
                   Kokkos::kokkos)
espresso_unit_test(SRC fene_test.cpp DEPENDS espresso::core Kokkos::kokkos)
espresso_unit_test(SRC Cell_boundary_test.cpp DEPENDS espresso::core
                   Kokkos::kokkos)
espresso_unit_test(SRC HaloPlan_test.cpp DEPENDS espresso::core Kokkos::kokkos
                   Boost::mpi MPI::MPI_CXX)
espresso_unit_test(SRC RegularHaloPlan_test.cpp DEPENDS espresso::core
                   Kokkos::kokkos Boost::mpi MPI::MPI_CXX NUM_PROC 4)
espresso_unit_test_executable(
  NAME HaloExchange_test SRC HaloExchange_test.cpp DEPENDS espresso::core
  Kokkos::kokkos Boost::mpi MPI::MPI_CXX)
foreach(TEST_NUM_PROC 1 2)
  if(${TEST_NUM_PROC} LESS_EQUAL ${ESPRESSO_TEST_NP})
    espresso_unit_test_register(
      NAME HaloExchange_test_${TEST_NUM_PROC}_mpi_ranks TARGET
      HaloExchange_test NUM_PROC ${TEST_NUM_PROC})
  endif()
endforeach()
# Register the validator test on both 4 ranks (exercises face/edge/corner peers
# on a non-trivial node grid) and 1 rank (where make_halo_plan is intentionally
# empty and the halo-layer cells are unreferenced ghosts).
espresso_unit_test_executable(
  NAME HaloPlanValidator_test SRC HaloPlanValidator_test.cpp DEPENDS
  espresso::core Kokkos::kokkos Boost::mpi MPI::MPI_CXX)
foreach(TEST_NUM_PROC 1 4)
  if(${TEST_NUM_PROC} LESS_EQUAL ${ESPRESSO_TEST_NP})
    espresso_unit_test_register(
      NAME HaloPlanValidator_test_${TEST_NUM_PROC}_mpi_ranks TARGET
      HaloPlanValidator_test NUM_PROC ${TEST_NUM_PROC})
  endif()
endforeach()
espresso_unit_test(SRC TabulatedPotential_test.cpp DEPENDS espresso::core
                   Kokkos::kokkos)
if(ESPRESSO_BUILD_WITH_WALBERLA)
  espresso_unit_test(SRC lb_particle_coupling_test.cpp DEPENDS espresso::core
                     Kokkos::kokkos Boost::mpi MPI::MPI_CXX NUM_PROC 2)
  espresso_unit_test(SRC ek_interface_test.cpp DEPENDS espresso::core
                     Kokkos::kokkos Boost::mpi MPI::MPI_CXX NUM_PROC 2)
  target_link_libraries(
    lb_particle_coupling_test PRIVATE espresso::walberla
                                      espresso::walberla::compiler_flags)
  target_link_libraries(
    ek_interface_test PRIVATE espresso::walberla
                              espresso::walberla::compiler_flags)
endif()

if(ESPRESSO_BUILD_WITH_FFTW)
  espresso_unit_test_executable(
    NAME p3m_test_matrix SRC p3m_test.cpp DEPENDS espresso::core Boost::mpi
    MPI::MPI_CXX
    $<$<BOOL:${ESPRESSO_BUILD_WITH_FFTW}>:espresso::p3m::compiler_flags>
    $<$<BOOL:${ESPRESSO_BUILD_WITH_FFTW}>:Heffte::Heffte> Cabana::Core
    Kokkos::kokkos)
  foreach(TEST_NUM_PROC 1 3 4)
    foreach(TEST_NUM_THREADS 1 2)
      if(${TEST_NUM_THREADS} LESS_EQUAL ${ESPRESSO_TEST_NT})
        if(${TEST_NUM_PROC} LESS_EQUAL ${ESPRESSO_TEST_NP})
          set(TEST_NAME
              p3m_test_${TEST_NUM_PROC}_mpi_ranks_${TEST_NUM_THREADS}_omp_threads
          )
          espresso_unit_test_register(
            NAME ${TEST_NAME} TARGET p3m_test_matrix NUM_PROC ${TEST_NUM_PROC}
            NUM_THREADS ${TEST_NUM_THREADS})
        endif()
      endif()
    endforeach()
  endforeach()
  espresso_unit_test(SRC fft_test.cpp DEPENDS espresso::utils espresso::core
                     Kokkos::kokkos Kokkos::kokkos)
  espresso_unit_test(SRC math_test.cpp DEPENDS espresso::utils espresso::core
                     Kokkos::kokkos)
  set_property(
    SOURCE math_test.cpp APPEND
    PROPERTY COMPILE_OPTIONS
             $<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:--diag-suppress=20208>)
endif()

if(ESPRESSO_BUILD_WITH_CUDA)
  espresso_unit_test(SRC cuda_test.cu DEPENDS espresso::utils espresso::core
                     Kokkos::kokkos CUDA::cudart)
endif()
set_property(
  SOURCE random_test.cpp thermostats_test.cpp
  APPEND
  PROPERTY
    COMPILE_OPTIONS
    $<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:--diag-suppress=20014>
    $<$<AND:$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>,$<STREQUAL:${CMAKE_CUDA_HOST_COMPILER_ID},NVHPC>>:--diag_suppress=declared_but_not_referenced>
)
set_property(
  SOURCE EspressoSystemStandAlone_test.cpp
  APPEND
  PROPERTY
    COMPILE_OPTIONS
    $<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:--diag-suppress=128>
    $<$<AND:$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>,$<STREQUAL:${CMAKE_CUDA_HOST_COMPILER_ID},NVHPC>>:--diag_suppress=loop_not_reachable>
)
set_property(
  SOURCE p3m_test.cpp
  APPEND
  PROPERTY
    COMPILE_OPTIONS
    $<$<AND:$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>,$<STREQUAL:${CMAKE_CUDA_HOST_COMPILER_ID},NVHPC>>:--diag_suppress=code_is_unreachable>
)
