#
# Copyright (C) 2010-2024 The ESPResSo project
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

espresso_unit_test(SRC Vector_test.cpp DEPENDS espresso::utils)
espresso_unit_test(SRC Factory_test.cpp DEPENDS espresso::utils)
espresso_unit_test(SRC NumeratedContainer_test.cpp DEPENDS espresso::utils)
espresso_unit_test(SRC keys_test.cpp DEPENDS espresso::utils)
espresso_unit_test(SRC Cache_test.cpp DEPENDS espresso::utils)
espresso_unit_test(SRC histogram_test.cpp DEPENDS espresso::utils)
espresso_unit_test(SRC accumulator_test.cpp DEPENDS espresso::utils
                   Boost::serialization)
espresso_unit_test(SRC int_pow_test.cpp DEPENDS espresso::utils)
espresso_unit_test(SRC sgn_test.cpp DEPENDS espresso::utils)
espresso_unit_test(SRC AS_erfc_part_test.cpp DEPENDS espresso::utils)
espresso_unit_test(SRC permute_ifield_test.cpp DEPENDS espresso::utils)
espresso_unit_test(SRC vec_rotate_test.cpp DEPENDS espresso::utils)
espresso_unit_test(SRC tensor_product_test.cpp DEPENDS espresso::utils)
espresso_unit_test(SRC linear_interpolation_test.cpp DEPENDS espresso::utils)
espresso_unit_test(SRC interpolation_gradient_test.cpp DEPENDS espresso::utils)
espresso_unit_test(SRC interpolation_test.cpp DEPENDS espresso::utils)
espresso_unit_test(SRC bspline_test.cpp DEPENDS espresso::utils)
espresso_unit_test(SRC matrix_vector_product_test.cpp DEPENDS espresso::utils)
espresso_unit_test(SRC index_test.cpp DEPENDS espresso::utils)
espresso_unit_test(SRC tuple_test.cpp DEPENDS espresso::utils)
espresso_unit_test(SRC Array_test.cpp DEPENDS Boost::serialization
                   espresso::utils)
espresso_unit_test(SRC contains_test.cpp DEPENDS espresso::utils)
espresso_unit_test(SRC Counter_test.cpp DEPENDS espresso::utils
                   Boost::serialization)
espresso_unit_test(SRC RunningAverage_test.cpp DEPENDS espresso::utils)
espresso_unit_test(SRC for_each_pair_test.cpp DEPENDS espresso::utils)
espresso_unit_test(SRC raster_test.cpp DEPENDS espresso::utils)
espresso_unit_test(SRC make_lin_space_test.cpp DEPENDS espresso::utils)
espresso_unit_test(SRC sampling_test.cpp DEPENDS espresso::utils)
espresso_unit_test(SRC coordinate_transformation_test.cpp DEPENDS
                   espresso::utils)
espresso_unit_test(SRC cylindrical_transformation_test.cpp DEPENDS
                   espresso::utils)
espresso_unit_test(SRC rotation_matrix_test.cpp DEPENDS espresso::utils)
espresso_unit_test(SRC quaternion_test.cpp DEPENDS espresso::utils)
espresso_unit_test(SRC mask_test.cpp DEPENDS espresso::utils)
espresso_unit_test(SRC type_traits_test.cpp DEPENDS espresso::utils)
espresso_unit_test(SRC uniform_test.cpp DEPENDS espresso::utils)
espresso_unit_test(SRC memcpy_archive_test.cpp DEPENDS espresso::utils)
espresso_unit_test(SRC triangle_functions_test.cpp DEPENDS espresso::utils)
espresso_unit_test(SRC Bag_test.cpp DEPENDS espresso::utils
                   Boost::serialization)
espresso_unit_test(SRC integral_parameter_test.cpp DEPENDS espresso::utils)
espresso_unit_test(SRC flatten_test.cpp DEPENDS espresso::utils)
espresso_unit_test(SRC pack_test.cpp DEPENDS Boost::serialization
                   espresso::utils)
espresso_unit_test(SRC unordered_map_test.cpp DEPENDS Boost::serialization
                   espresso::utils)
espresso_unit_test(SRC u32_to_u64_test.cpp DEPENDS espresso::utils NUM_PROC 1)
espresso_unit_test(SRC gather_buffer_test.cpp DEPENDS espresso::utils::mpi
                   Boost::mpi MPI::MPI_CXX NUM_PROC 4)
espresso_unit_test(SRC scatter_buffer_test.cpp DEPENDS espresso::utils::mpi
                   Boost::mpi MPI::MPI_CXX NUM_PROC 4)
espresso_unit_test(SRC all_compare_test.cpp DEPENDS espresso::utils::mpi
                   Boost::mpi MPI::MPI_CXX NUM_PROC 3)
espresso_unit_test(SRC gatherv_test.cpp DEPENDS espresso::utils::mpi Boost::mpi
                   MPI::MPI_CXX NUM_PROC 3)
espresso_unit_test(SRC iall_gatherv_test.cpp DEPENDS espresso::utils::mpi
                   Boost::mpi MPI::MPI_CXX NUM_PROC 3)
espresso_unit_test(SRC sendrecv_test.cpp DEPENDS espresso::utils::mpi
                   Boost::mpi MPI::MPI_CXX espresso::utils NUM_PROC 3)
espresso_unit_test(SRC serialization_test.cpp DEPENDS espresso::utils
                   Boost::serialization Boost::mpi MPI::MPI_CXX NUM_PROC 1)
espresso_unit_test(SRC matrix_test.cpp DEPENDS espresso::utils
                   Boost::serialization NUM_PROC 1)
espresso_unit_test(SRC orthonormal_vec_test.cpp DEPENDS espresso::utils
                   Boost::serialization NUM_PROC 1)
espresso_unit_test(SRC reduce_optional_test.cpp DEPENDS espresso::utils::mpi
                   Boost::mpi MPI::MPI_CXX NUM_PROC 4)
