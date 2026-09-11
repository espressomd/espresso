#
# Copyright (C) 2009-2026 The ESPResSo project
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
  bond_error.cpp
  bonds.cpp
  cells.cpp
  communication.cpp
  dpd.cpp
  energy.cpp
  errorhandling.cpp
  forces.cpp
  forces_init.cpp
  ghosts/HaloExchange.cpp
  ghosts/HaloPlanValidator.cpp
  ghosts/particle_packing.cpp
  integrate.cpp
  npt.cpp
  particle_node.cpp
  particle_store/MigrationPack.cpp
  particle_store/ParticleStore.cpp
  polymer.cpp
  pressure.cpp
  propagation.cpp
  rattle.cpp
  rotation.cpp
  short_range_verlet.cpp
  Observable_stat.cpp
  thermostat.cpp
  tuning.cpp
  virtual_sites.cpp
  exclusions.cpp
  system/GpuParticleData.cpp
  system/System.cpp
  PartCfg.cpp
  ExclusionRadius.cpp
  TabulatedPotential.cpp)
add_library(espresso::core ALIAS espresso_core)
espresso_set_common_target_properties(espresso_core)

if(ESPRESSO_BUILD_WITH_CUDA)
  target_sources(
    espresso_core
    PRIVATE cuda/init.cpp cuda/common_cuda.cu cuda/init_cuda.cu
            cuda/CudaHostAllocator.cu forces_lj_device.cu
            magnetostatics/dipolar_direct_sum_gpu_cuda.cu
            electrostatics/p3m_gpu_cuda.cu electrostatics/p3m_gpu_error_cuda.cu
            system/GpuParticleData_cuda.cu)
  target_link_libraries(espresso_core PRIVATE CUDA::cuda_driver CUDA::cudart
                                              CUDA::cufft)
  target_include_directories(
    espresso_core
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
            ${CMAKE_CURRENT_SOURCE_DIR}/electrostatics
            ${CMAKE_CURRENT_SOURCE_DIR}/magnetostatics)
endif()

if(ESPRESSO_BUILD_WITH_NLOPT)
  target_include_directories(espresso_core PRIVATE "${nlopt_SOURCE_DIR}/src")
endif()

espresso_set_kokkos_properties_on_source_files(espresso_core)
install(TARGETS espresso_core
        LIBRARY DESTINATION ${ESPRESSO_INSTALL_PYTHON}/espressomd)

target_link_libraries(
  espresso_core
  PUBLIC Cabana::Core Kokkos::kokkos OpenMP::OpenMP_CXX
         $<$<BOOL:${ESPRESSO_BUILD_WITH_CUDA}>:OpenMP::OpenMP_CUDA>
  PRIVATE espresso::config espresso::utils::mpi espresso::shapes
          espresso::compiler_flags
          $<$<BOOL:${ESPRESSO_BUILD_WITH_WALBERLA}>:espresso::walberla>
          $<$<BOOL:${ESPRESSO_BUILD_WITH_NLOPT}>:nlopt>
          $<$<BOOL:${ESPRESSO_BUILD_WITH_GSL}>:GSL::gsl>
  PUBLIC espresso::utils MPI::MPI_CXX Random123 espresso::particle_observables
         Boost::serialization Boost::mpi espresso::instrumentation)

target_include_directories(espresso_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

if(ESPRESSO_BUILD_WITH_FFTW)
  add_subdirectory(p3m)
  add_subdirectory(fft)
endif()

add_subdirectory(accumulators)
add_subdirectory(analysis)
add_subdirectory(bond_breakage)
add_subdirectory(bonded_interactions)
add_subdirectory(cell_system)
add_subdirectory(cluster_analysis)
add_subdirectory(collision_detection)
add_subdirectory(constraints)
add_subdirectory(ek)
add_subdirectory(electrostatics)
add_subdirectory(error_handling)
add_subdirectory(galilei)
add_subdirectory(immersed_boundary)
add_subdirectory(integrators)
add_subdirectory(io)
add_subdirectory(lb)
add_subdirectory(magnetostatics)
add_subdirectory(nonbonded_interactions)
add_subdirectory(object-in-fluid)
add_subdirectory(observables)
add_subdirectory(scafacos)
add_subdirectory(stokesian_dynamics)
add_subdirectory(virtual_sites)

if(ESPRESSO_BUILD_TESTS)
  add_subdirectory(unit_tests)
endif()

if(ESPRESSO_BUILD_WITH_HDF5 AND ESPRESSO_BUILD_WITH_COVERAGE
   AND NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  set_property(SOURCE io/writer/h5md_core.cpp APPEND
               PROPERTY COMPILE_OPTIONS -felide-constructors)
endif()
# silence inconsistent host/device qualifiers in Boost accumulators and combine
set_property(
  SOURCE electrostatics/coulomb.cpp observables/RDF.cpp APPEND
  PROPERTY COMPILE_OPTIONS
           $<$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>:--diag-suppress=20014>)
set_property(
  SOURCE electrostatics/coulomb.cpp
  APPEND
  PROPERTY
    COMPILE_OPTIONS
    $<$<AND:$<COMPILE_LANG_AND_ID:CUDA,NVIDIA>,$<STREQUAL:${CMAKE_CUDA_HOST_COMPILER_ID},NVHPC>>:--diag_suppress=code_is_unreachable>
)

espresso_set_kokkos_properties_on_source_files(espresso_core)
if(ESPRESSO_BUILD_WITH_FFTW)
  espresso_set_kokkos_properties_on_source_files(espresso_p3m)
endif()
