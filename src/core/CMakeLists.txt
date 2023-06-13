#
# Copyright (C) 2009-2022 The ESPResSo project
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

add_library(
  espresso_core SHARED
  accumulators.cpp
  bond_error.cpp
  cells.cpp
  collision.cpp
  communication.cpp
  constraints.cpp
  dpd.cpp
  energy.cpp
  errorhandling.cpp
  EspressoSystemInterface.cpp
  forcecap.cpp
  forces.cpp
  ghosts.cpp
  grid.cpp
  immersed_boundaries.cpp
  interactions.cpp
  event.cpp
  integrate.cpp
  npt.cpp
  partCfg_global.cpp
  particle_node.cpp
  polymer.cpp
  pressure.cpp
  rattle.cpp
  rotate_system.cpp
  rotation.cpp
  Observable_stat.cpp
  SystemInterface.cpp
  thermostat.cpp
  tuning.cpp
  virtual_sites.cpp
  exclusions.cpp
  PartCfg.cpp
  EspressoSystemStandAlone.cpp
  TabulatedPotential.cpp)
add_library(espresso::core ALIAS espresso_core)
set_target_properties(espresso_core PROPERTIES CXX_CLANG_TIDY
                                               "${ESPRESSO_CXX_CLANG_TIDY}")

if(ESPRESSO_BUILD_WITH_CUDA)
  target_sources(espresso_core PRIVATE cuda_init.cpp cuda_interface.cpp)
  espresso_add_gpu_library(
    espresso_cuda SHARED cuda_common_cuda.cu cuda_init_cuda.cu
    CudaHostAllocator.cu magnetostatics/barnes_hut_gpu_cuda.cu
    magnetostatics/dipolar_direct_sum_gpu_cuda.cu
    electrostatics/mmm1d_gpu_cuda.cu electrostatics/p3m_gpu_cuda.cu
    electrostatics/p3m_gpu_error_cuda.cu EspressoSystemInterface_cuda.cu)
  add_library(espresso::cuda ALIAS espresso_cuda)
  target_link_libraries(
    espresso_cuda PRIVATE CUDA::cuda_driver CUDA::cudart CUDA::cufft
                          espresso::config espresso::utils espresso::shapes)
  target_include_directories(
    espresso_cuda
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
            ${CMAKE_CURRENT_SOURCE_DIR}/electrostatics
            ${CMAKE_CURRENT_SOURCE_DIR}/magnetostatics)
  set_target_properties(espresso_cuda PROPERTIES CUDA_CLANG_TIDY
                                                 "${ESPRESSO_CUDA_CLANG_TIDY}")
  target_link_libraries(espresso_core PUBLIC espresso::cuda)
  install(TARGETS espresso_cuda
          LIBRARY DESTINATION ${ESPRESSO_INSTALL_PYTHON}/espressomd)
endif()

install(TARGETS espresso_core
        LIBRARY DESTINATION ${ESPRESSO_INSTALL_PYTHON}/espressomd)

target_link_libraries(
  espresso_core PRIVATE espresso::config espresso::utils::mpi espresso::shapes
                        espresso::profiler espresso::cpp_flags
  PUBLIC espresso::utils MPI::MPI_CXX Random123 espresso::particle_observables
         Boost::serialization Boost::mpi)

target_include_directories(espresso_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_subdirectory(accumulators)
add_subdirectory(analysis)
add_subdirectory(bond_breakage)
add_subdirectory(bonded_interactions)
add_subdirectory(cell_system)
add_subdirectory(cluster_analysis)
add_subdirectory(constraints)
add_subdirectory(electrostatics)
add_subdirectory(error_handling)
add_subdirectory(galilei)
add_subdirectory(grid_based_algorithms)
add_subdirectory(immersed_boundary)
add_subdirectory(integrators)
add_subdirectory(io)
add_subdirectory(magnetostatics)
add_subdirectory(nonbonded_interactions)
add_subdirectory(object-in-fluid)
add_subdirectory(observables)
add_subdirectory(p3m)
add_subdirectory(reaction_methods)
add_subdirectory(scafacos)
add_subdirectory(stokesian_dynamics)
add_subdirectory(virtual_sites)

if(ESPRESSO_BUILD_TESTS)
  add_subdirectory(unit_tests)
endif()
