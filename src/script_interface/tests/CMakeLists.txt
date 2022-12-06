#
# Copyright (C) 2020-2022 The ESPResSo project
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

unit_test(NAME ObjectHandle_test SRC ObjectHandle_test.cpp DEPENDS
          espresso::script_interface)
unit_test(NAME AutoParameters_test SRC AutoParameters_test.cpp DEPENDS
          espresso::script_interface)
unit_test(NAME AutoParameter_test SRC AutoParameter_test.cpp DEPENDS
          espresso::script_interface)
unit_test(NAME Variant_test SRC Variant_test.cpp DEPENDS
          espresso::script_interface)
unit_test(NAME get_value_test SRC get_value_test.cpp DEPENDS
          espresso::script_interface)
unit_test(NAME None_test SRC None_test.cpp DEPENDS espresso::script_interface)
unit_test(NAME reduction_test SRC reduction_test.cpp DEPENDS
          espresso::script_interface Boost::mpi MPI::MPI_CXX NUM_PROC 4)
unit_test(NAME LocalContext_test SRC LocalContext_test.cpp DEPENDS
          espresso::script_interface Boost::mpi MPI::MPI_CXX NUM_PROC 1)
unit_test(NAME GlobalContext_test SRC GlobalContext_test.cpp DEPENDS
          espresso::script_interface Boost::mpi MPI::MPI_CXX NUM_PROC 2)
unit_test(NAME Exception_test SRC Exception_test.cpp DEPENDS
          espresso::script_interface)
unit_test(NAME ParallelExceptionHandler_test SRC
          ParallelExceptionHandler_test.cpp DEPENDS espresso::script_interface
          espresso::core Boost::mpi MPI::MPI_CXX NUM_PROC 2)
unit_test(NAME packed_variant_test SRC packed_variant_test.cpp DEPENDS
          espresso::script_interface)
unit_test(NAME ObjectList_test SRC ObjectList_test.cpp DEPENDS
          espresso::script_interface espresso::core Boost::mpi)
unit_test(NAME ObjectMap_test SRC ObjectMap_test.cpp DEPENDS
          espresso::script_interface espresso::core Boost::mpi)
unit_test(NAME serialization_mpi_guard_test SRC
          serialization_mpi_guard_test.cpp DEPENDS espresso::script_interface
          Boost::mpi MPI::MPI_CXX NUM_PROC 2)
unit_test(NAME Accumulators_test SRC Accumulators_test.cpp DEPENDS
          espresso::script_interface espresso::core)
unit_test(NAME Constraints_test SRC Constraints_test.cpp DEPENDS
          espresso::script_interface espresso::core)
unit_test(NAME Actors_test SRC Actors_test.cpp DEPENDS
          espresso::script_interface espresso::core)
