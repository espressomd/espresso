# Copyright (C) 2010-2019 The ESPResSo project
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

include(unit_test)

# Add tests here

unit_test(NAME RuntimeError_test SRC RuntimeError_test.cpp DEPENDS
          Boost::serialization)
unit_test(NAME RuntimeErrorCollector_test SRC RuntimeErrorCollector_test.cpp
          DEPENDS EspressoCore Boost::mpi MPI::MPI_CXX)

unit_test(NAME ScriptInterface_test SRC ScriptInterface_test.cpp DEPENDS
          ScriptInterface)
unit_test(NAME MpiCallbacks_test SRC MpiCallbacks_test.cpp DEPENDS utils
          Boost::mpi MPI::MPI_CXX NUM_PROC 2)
unit_test(NAME ParallelScriptInterface_test SRC
          ParallelScriptInterface_test.cpp DEPENDS ScriptInterface Boost::mpi
          MPI::MPI_CXX NUM_PROC 2)
unit_test(NAME AutoParameters_test SRC AutoParameters_test.cpp DEPENDS
          ScriptInterface)
unit_test(NAME AutoParameter_test SRC AutoParameter_test.cpp DEPENDS
          ScriptInterface)
unit_test(NAME Variant_test SRC Variant_test.cpp DEPENDS ScriptInterface)
unit_test(NAME ParticleIterator_test SRC ParticleIterator_test.cpp DEPENDS
          utils)
unit_test(NAME link_cell_test SRC link_cell_test.cpp DEPENDS utils)
unit_test(NAME verlet_ia_test SRC verlet_ia_test.cpp DEPENDS utils)
unit_test(NAME Particle_test SRC Particle_test.cpp DEPENDS utils
          Boost::serialization)
unit_test(NAME get_value SRC get_value_test.cpp DEPENDS ScriptInterface)
unit_test(NAME field_coupling_couplings SRC field_coupling_couplings_test.cpp
          DEPENDS utils)
unit_test(NAME field_coupling_fields SRC field_coupling_fields_test.cpp DEPENDS
          utils)
unit_test(NAME field_coupling_force_field SRC
          field_coupling_force_field_test.cpp DEPENDS utils)
unit_test(NAME periodic_fold_test SRC periodic_fold_test.cpp)
unit_test(NAME None_test SRC None_test.cpp DEPENDS ScriptInterface)
unit_test(NAME grid_test SRC grid_test.cpp DEPENDS EspressoCore)
unit_test(NAME BoxGeometry_test SRC BoxGeometry_test.cpp DEPENDS EspressoCore)
unit_test(NAME LocalBox_test SRC LocalBox_test.cpp DEPENDS EspressoCore)
unit_test(NAME thermostats_test SRC thermostats_test.cpp DEPENDS EspressoCore)
unit_test(NAME random_test SRC random_test.cpp DEPENDS utils Random123)
unit_test(NAME BondList_test SRC BondList_test.cpp DEPENDS EspressoCore)
