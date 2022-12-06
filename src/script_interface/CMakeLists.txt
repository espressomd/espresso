#
# Copyright (C) 2016-2022 The ESPResSo project
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
  espresso_script_interface SHARED
  initialize.cpp ObjectHandle.cpp object_container_mpi_guard.cpp
  GlobalContext.cpp ContextManager.cpp ParallelExceptionHandler.cpp)
add_library(espresso::script_interface ALIAS espresso_script_interface)

add_subdirectory(accumulators)
add_subdirectory(analysis)
add_subdirectory(bond_breakage)
add_subdirectory(cell_system)
add_subdirectory(cluster_analysis)
add_subdirectory(code_info)
add_subdirectory(collision_detection)
add_subdirectory(constraints)
add_subdirectory(electrostatics)
add_subdirectory(galilei)
add_subdirectory(h5md)
add_subdirectory(integrators)
add_subdirectory(interactions)
add_subdirectory(lees_edwards)
add_subdirectory(magnetostatics)
add_subdirectory(math)
add_subdirectory(mpiio)
add_subdirectory(observables)
add_subdirectory(pair_criteria)
add_subdirectory(particle_data)
add_subdirectory(reaction_methods)
add_subdirectory(scafacos)
add_subdirectory(shapes)
add_subdirectory(system)
add_subdirectory(virtual_sites)
add_subdirectory(walberla)

install(TARGETS espresso_script_interface
        LIBRARY DESTINATION ${ESPRESSO_INSTALL_PYTHON}/espressomd)

target_link_libraries(
  espresso_script_interface PRIVATE espresso::config espresso::core
  PUBLIC espresso::utils MPI::MPI_CXX Boost::mpi espresso::shapes
  PRIVATE espresso::cpp_flags)

target_include_directories(espresso_script_interface
                           PUBLIC ${CMAKE_SOURCE_DIR}/src)

if(ESPRESSO_BUILD_TESTS)
  add_subdirectory(tests)
endif()
