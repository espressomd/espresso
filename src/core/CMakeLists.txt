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

set(Espresso_core_SRC
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
    particle_data.cpp
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
    EspressoSystemStandAlone.cpp)

set(Espresso_cuda_SRC
    cuda_common_cuda.cu
    cuda_init.cpp
    cuda_init_cuda.cu
    cuda_interface.cpp
    CudaHostAllocator.cu
    magnetostatics/barnes_hut_gpu_cuda.cu
    magnetostatics/dds_gpu_cuda.cu
    electrostatics/mmm1d_gpu_cuda.cu
    electrostatics/p3m_gpu_cuda.cu
    electrostatics/p3m_gpu_error_cuda.cu
    EspressoSystemInterface_cuda.cu
    grid_based_algorithms/electrokinetics_cuda.cu
    grid_based_algorithms/lbgpu_cuda.cu
    grid_based_algorithms/fd-electrostatics_cuda.cu
    grid_based_algorithms/electrokinetics.cpp
    grid_based_algorithms/lbgpu.cpp
    virtual_sites/lb_inertialess_tracers_cuda.cu)

if(CUDA)
  add_gpu_library(Espresso_core SHARED ${Espresso_core_SRC}
                  ${Espresso_cuda_SRC})
else()
  add_library(Espresso_core SHARED ${Espresso_core_SRC})
endif()
add_library(Espresso::core ALIAS Espresso_core)

install(TARGETS Espresso_core LIBRARY DESTINATION ${PYTHON_INSTDIR}/espressomd)

target_link_libraries(
  Espresso_core
  PRIVATE Espresso::config Espresso::shapes Espresso::profiler
          $<$<BOOL:${SCAFACOS}>:Espresso::scafacos> Espresso::cpp_flags
  PUBLIC Espresso::utils MPI::MPI_CXX Random123 Espresso::particle_observables
         Boost::serialization Boost::mpi "$<$<BOOL:${H5MD}>:${HDF5_LIBRARIES}>"
         $<$<BOOL:${H5MD}>:Boost::filesystem> $<$<BOOL:${H5MD}>:h5xx>
         $<$<BOOL:${FFTW3_FOUND}>:FFTW3::FFTW3>)

target_include_directories(
  Espresso_core
  PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
         "$<$<BOOL:${H5MD}>:${CMAKE_CURRRENT_SOURCE_DIR}/io/writer>"
         "$<$<BOOL:${H5MD}>:${HDF5_INCLUDE_DIRS}>")

target_compile_definitions(Espresso_core PUBLIC $<$<BOOL:${H5MD}>:H5XX_USE_MPI>)

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
add_subdirectory(virtual_sites)

if(WITH_TESTS)
  add_subdirectory(unit_tests)
endif(WITH_TESTS)

if(STOKESIAN_DYNAMICS)
  add_subdirectory(stokesian_dynamics)
  target_link_libraries(Espresso_core PRIVATE StokesianDynamics::sd_cpu)
endif()
