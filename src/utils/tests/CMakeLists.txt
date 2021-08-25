include(unit_test)

unit_test(NAME abs_test SRC abs_test.cpp DEPENDS EspressoUtils)
unit_test(NAME Vector_test SRC Vector_test.cpp DEPENDS EspressoUtils)
unit_test(NAME Factory_test SRC Factory_test.cpp DEPENDS EspressoUtils)
unit_test(NAME NumeratedContainer_test SRC NumeratedContainer_test.cpp DEPENDS
          EspressoUtils)
unit_test(NAME keys_test SRC keys_test.cpp DEPENDS EspressoUtils)
unit_test(NAME Cache_test SRC Cache_test.cpp DEPENDS EspressoUtils)
unit_test(NAME histogram SRC histogram.cpp DEPENDS EspressoUtils)
unit_test(NAME accumulator SRC accumulator.cpp DEPENDS EspressoUtils
          Boost::serialization)
unit_test(NAME int_pow SRC int_pow_test.cpp DEPENDS EspressoUtils)
unit_test(NAME sgn SRC sgn_test.cpp DEPENDS EspressoUtils)
unit_test(NAME AS_erfc_part SRC AS_erfc_part_test.cpp DEPENDS EspressoUtils)
unit_test(NAME sinc SRC sinc_test.cpp DEPENDS EspressoUtils)
unit_test(NAME as_const SRC as_const_test.cpp DEPENDS EspressoUtils)
unit_test(NAME permute_ifield_test SRC permute_ifield_test.cpp DEPENDS
          EspressoUtils)
unit_test(NAME vec_rotate SRC vec_rotate_test.cpp DEPENDS EspressoUtils)
unit_test(NAME tensor_product SRC tensor_product_test.cpp DEPENDS EspressoUtils)
unit_test(NAME linear_interpolation SRC linear_interpolation_test.cpp DEPENDS
          EspressoUtils)
unit_test(NAME interpolation_gradient SRC interpolation_gradient_test.cpp
          DEPENDS EspressoUtils)
unit_test(NAME interpolation SRC interpolation_test.cpp DEPENDS EspressoUtils)
unit_test(NAME bspline_test SRC bspline_test.cpp DEPENDS EspressoUtils)
unit_test(NAME Span_test SRC Span_test.cpp DEPENDS EspressoUtils)
unit_test(NAME matrix_vector_product SRC matrix_vector_product.cpp DEPENDS
          EspressoUtils)
unit_test(NAME index_test SRC index_test.cpp DEPENDS EspressoUtils)
unit_test(NAME tuple_test SRC tuple_test.cpp DEPENDS EspressoUtils)
unit_test(NAME Array_test SRC Array_test.cpp DEPENDS Boost::serialization
          EspressoUtils)
unit_test(NAME contains_test SRC contains_test.cpp DEPENDS EspressoUtils)
unit_test(NAME Counter_test SRC Counter_test.cpp DEPENDS EspressoUtils
          Boost::serialization)
unit_test(NAME RunningAverage_test SRC RunningAverage_test.cpp DEPENDS
          EspressoUtils)
unit_test(NAME for_each_pair_test SRC for_each_pair_test.cpp DEPENDS
          EspressoUtils)
unit_test(NAME raster_test SRC raster_test.cpp DEPENDS EspressoUtils)
unit_test(NAME make_lin_space_test SRC make_lin_space_test.cpp DEPENDS
          EspressoUtils)
unit_test(NAME sampling_test SRC sampling_test.cpp DEPENDS EspressoUtils)
unit_test(NAME coordinate_transformation_test SRC coordinate_transformation.cpp
          DEPENDS EspressoUtils)
unit_test(NAME cylindrical_transformation_test SRC
          cylindrical_transformation.cpp DEPENDS EspressoUtils)
unit_test(NAME rotation_matrix_test SRC rotation_matrix_test.cpp DEPENDS
          EspressoUtils)
unit_test(NAME quaternion_test SRC quaternion_test.cpp DEPENDS EspressoUtils)
unit_test(NAME mask_test SRC mask_test.cpp DEPENDS EspressoUtils)
unit_test(NAME type_traits_test SRC type_traits_test.cpp DEPENDS EspressoUtils)
unit_test(NAME uniform_test SRC uniform_test.cpp DEPENDS EspressoUtils)
unit_test(NAME memcpy_archive_test SRC memcpy_archive_test.cpp DEPENDS
          EspressoUtils)
unit_test(NAME triangle_functions_test SRC triangle_functions_test.cpp DEPENDS
          EspressoUtils)
unit_test(NAME Bag_test SRC Bag_test.cpp DEPENDS EspressoUtils
          Boost::serialization)
unit_test(NAME integral_parameter_test SRC integral_parameter_test.cpp DEPENDS
          EspressoUtils)
unit_test(NAME flatten_test SRC flatten_test.cpp DEPENDS EspressoUtils)
unit_test(NAME pack_test SRC pack_test.cpp DEPENDS Boost::serialization
          EspressoUtils)
unit_test(NAME unordered_map_test SRC unordered_map_test.cpp DEPENDS
          Boost::serialization EspressoUtils)
unit_test(NAME u32_to_u64_test SRC u32_to_u64_test.cpp DEPENDS EspressoUtils
          NUM_PROC 1)
unit_test(NAME gather_buffer_test SRC gather_buffer_test.cpp DEPENDS
          EspressoUtils Boost::mpi MPI::MPI_CXX NUM_PROC 4)
unit_test(NAME scatter_buffer_test SRC scatter_buffer_test.cpp DEPENDS
          EspressoUtils Boost::mpi MPI::MPI_CXX NUM_PROC 4)
unit_test(NAME all_compare_test SRC all_compare_test.cpp DEPENDS EspressoUtils
          Boost::mpi MPI::MPI_CXX NUM_PROC 3)
unit_test(NAME gatherv_test SRC gatherv_test.cpp DEPENDS EspressoUtils
          Boost::mpi MPI::MPI_CXX NUM_PROC 3)
unit_test(NAME sendrecv_test SRC sendrecv_test.cpp DEPENDS EspressoUtils
          Boost::mpi MPI::MPI_CXX EspressoUtils NUM_PROC 3)
unit_test(NAME matrix_test SRC matrix_test.cpp DEPENDS EspressoUtils
          Boost::serialization NUM_PROC 1)
unit_test(NAME orthonormal_vec_test SRC orthonormal_vec_test.cpp DEPENDS
          EspressoUtils Boost::serialization NUM_PROC 1)
