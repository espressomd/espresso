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
  Espresso_script_interface SHARED
  initialize.cpp ObjectHandle.cpp object_container_mpi_guard.cpp
  GlobalContext.cpp ContextManager.cpp ParallelExceptionHandler.cpp)
add_library(Espresso::script_interface ALIAS Espresso_script_interface)

add_subdirectory(accumulators)
add_subdirectory(bond_breakage)
add_subdirectory(cell_system)
add_subdirectory(collision_detection)
add_subdirectory(constraints)
add_subdirectory(cluster_analysis)
add_subdirectory(electrostatics)
add_subdirectory(interactions)
add_subdirectory(lbboundaries)
add_subdirectory(lees_edwards)
add_subdirectory(magnetostatics)
add_subdirectory(virtual_sites)
add_subdirectory(observables)
add_subdirectory(pair_criteria)
add_subdirectory(mpiio)
add_subdirectory(scafacos)
add_subdirectory(shapes)
add_subdirectory(h5md)
add_subdirectory(reaction_methods)

install(TARGETS Espresso_script_interface
        LIBRARY DESTINATION ${PYTHON_INSTDIR}/espressomd)

target_link_libraries(
  Espresso_script_interface PRIVATE Espresso::config Espresso::core
  PUBLIC Espresso::utils MPI::MPI_CXX Boost::mpi Espresso::shapes
  PRIVATE Espresso::cpp_flags)

target_include_directories(Espresso_script_interface
                           PUBLIC ${CMAKE_SOURCE_DIR}/src)

add_subdirectory(tests)
