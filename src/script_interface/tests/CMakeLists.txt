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

espresso_unit_test(SRC ObjectHandle_test.cpp DEPENDS espresso::script_interface)
espresso_unit_test(SRC AutoParameters_test.cpp DEPENDS
                   espresso::script_interface)
espresso_unit_test(SRC AutoParameter_test.cpp DEPENDS
                   espresso::script_interface)
espresso_unit_test(SRC Variant_test.cpp DEPENDS espresso::script_interface)
espresso_unit_test(SRC get_value_test.cpp DEPENDS espresso::script_interface)
espresso_unit_test(SRC None_test.cpp DEPENDS espresso::script_interface)
espresso_unit_test(SRC reduction_test.cpp DEPENDS espresso::script_interface
                   Boost::mpi MPI::MPI_CXX NUM_PROC 4)
espresso_unit_test(SRC LocalContext_test.cpp DEPENDS espresso::script_interface
                   Boost::mpi MPI::MPI_CXX NUM_PROC 1)
espresso_unit_test(
  SRC GlobalContext_test.cpp DEPENDS espresso::script_interface Boost::mpi
  MPI::MPI_CXX NUM_PROC 2)
espresso_unit_test(SRC Exception_test.cpp DEPENDS espresso::script_interface)
espresso_unit_test(
  SRC ParallelExceptionHandler_test.cpp DEPENDS espresso::script_interface
  espresso::core Boost::mpi MPI::MPI_CXX NUM_PROC 2)
espresso_unit_test(SRC packed_variant_test.cpp DEPENDS
                   espresso::script_interface)
espresso_unit_test(SRC ObjectList_test.cpp DEPENDS espresso::script_interface
                   espresso::core Boost::mpi)
espresso_unit_test(SRC ObjectMap_test.cpp DEPENDS espresso::script_interface
                   espresso::core Boost::mpi)
espresso_unit_test(SRC Accumulators_test.cpp DEPENDS espresso::script_interface
                   espresso::core Boost::mpi MPI::MPI_CXX NUM_PROC 2)
espresso_unit_test(SRC Constraints_test.cpp DEPENDS espresso::script_interface
                   espresso::core)
espresso_unit_test(SRC Actors_test.cpp DEPENDS espresso::script_interface
                   espresso::core)
espresso_unit_test(
  SRC ConstantpHEnsemble_test.cpp DEPENDS espresso::core
  espresso::script_interface Boost::mpi MPI::MPI_CXX NUM_PROC 2)
espresso_unit_test(
  SRC ReactionEnsemble_test.cpp DEPENDS espresso::core
  espresso::script_interface Boost::mpi MPI::MPI_CXX NUM_PROC 2)
