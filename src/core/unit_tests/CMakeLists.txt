# Copyright (C) 2010-2018 The ESPResSo project
# Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
#   Max-Planck-Institute for Polymer Research, Theory Group

# This file is part of ESPResSo.

# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>. 

# Add tests here
unit_test(NAME Vector_test SRC Vector_test.cpp)
unit_test(NAME RuntimeError_test SRC RuntimeError_test.cpp DEPENDS Boost::serialization)
unit_test(NAME RunningAverage_test SRC RunningAverage_test.cpp)
unit_test(NAME RuntimeErrorCollector_test SRC RuntimeErrorCollector_test.cpp DEPENDS EspressoCore Boost::mpi MPI::MPI_CXX)

unit_test(NAME ScriptInterface_test SRC ScriptInterface_test.cpp DEPENDS EspressoScriptInterface)
unit_test(NAME Wall_test SRC Wall_test.cpp ../shapes/Wall.cpp)
unit_test(NAME Ellipsoid_test SRC Ellipsoid_test.cpp DEPENDS Shapes EspressoConfig)
unit_test(NAME Factory_test SRC Factory_test.cpp)
unit_test(NAME NumeratedContainer_test SRC NumeratedContainer_test.cpp)
unit_test(NAME MpiCallbacks_test SRC MpiCallbacks_test.cpp ../MpiCallbacks.cpp DEPENDS Boost::mpi MPI::MPI_CXX NUM_PROC 2)
unit_test(NAME ParallelScriptInterface_test SRC ParallelScriptInterface_test.cpp DEPENDS EspressoScriptInterface Boost::mpi MPI::MPI_CXX NUM_PROC 2)
unit_test(NAME gather_buffer_test SRC gather_buffer_test.cpp DEPENDS Boost::mpi MPI::MPI_CXX NUM_PROC 4)
unit_test(NAME scatter_buffer_test SRC scatter_buffer_test.cpp DEPENDS Boost::mpi MPI::MPI_CXX NUM_PROC 4)
unit_test(NAME AutoParameters_test SRC AutoParameters_test.cpp DEPENDS EspressoScriptInterface)
unit_test(NAME make_function_test SRC make_function_test.cpp)
unit_test(NAME AutoParameter_test SRC AutoParameter_test.cpp DEPENDS EspressoScriptInterface)
unit_test(NAME Variant_test SRC Variant_test.cpp DEPENDS EspressoScriptInterface)
unit_test(NAME ParticleIterator_test SRC ParticleIterator_test.cpp)
unit_test(NAME link_cell_test SRC link_cell_test.cpp)
unit_test(NAME verlet_ia_test SRC verlet_ia_test.cpp)
unit_test(NAME ParticleCache_test SRC ParticleCache_test.cpp ../MpiCallbacks.cpp DEPENDS Boost::mpi MPI::MPI_CXX Boost::serialization NUM_PROC 2)
unit_test(NAME SkipIterator_test SRC SkipIterator_test.cpp)
unit_test(NAME List_test SRC List_test.cpp DEPENDS Boost::serialization)
unit_test(NAME Particle_test SRC Particle_test.cpp DEPENDS Boost::serialization)
unit_test(NAME for_each_pair_test SRC for_each_pair_test.cpp)
unit_test(NAME all_compare_test SRC all_compare_test.cpp DEPENDS Boost::mpi MPI::MPI_CXX NUM_PROC 3)
unit_test(NAME None_test SRC None_test.cpp DEPENDS EspressoScriptInterface)
unit_test(NAME keys_test SRC keys_test.cpp)
unit_test(NAME Cache_test SRC Cache_test.cpp)
unit_test(NAME histogram SRC histogram.cpp)
unit_test(NAME get_value SRC get_value_test.cpp DEPENDS EspressoScriptInterface)
unit_test(NAME accumulator SRC accumulator.cpp)
unit_test(NAME strcat_alloc SRC strcat_alloc_test.cpp)
unit_test(NAME int_pow SRC int_pow_test.cpp)
unit_test(NAME sgn SRC sgn_test.cpp)
unit_test(NAME AS_erfc_part SRC AS_erfc_part_test.cpp)
unit_test(NAME sinc SRC sinc_test.cpp)
unit_test(NAME as_const SRC as_const_test.cpp)
unit_test(NAME permute_ifield_test SRC permute_ifield_test.cpp)
unit_test(NAME vec_rotate SRC vec_rotate_test.cpp)
unit_test(NAME tensor_product SRC tensor_product_test.cpp)
unit_test(NAME interpolation_gradient SRC interpolation_gradient_test.cpp)
unit_test(NAME interpolation SRC interpolation_test.cpp)
unit_test(NAME field_coupling_couplings SRC field_coupling_couplings_test.cpp)
unit_test(NAME field_coupling_fields SRC field_coupling_fields_test.cpp)
unit_test(NAME field_coupling_force_field SRC field_coupling_force_field_test.cpp)
unit_test(NAME Span_test SRC Span_test.cpp)
unit_test(NAME periodic_fold_test SRC periodic_fold_test.cpp)
unit_test(NAME sendrecv_test SRC sendrecv_test.cpp DEPENDS Boost::mpi MPI::MPI_CXX)
unit_test(NAME LbWalberla_test SRC LbWalberla.cpp DEPENDS EspressoCore Boost::mpi MPI::MPI_CXX MPI::MPI_C)
