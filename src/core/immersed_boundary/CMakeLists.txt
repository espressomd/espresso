file(GLOB ImmersedBoundary_SRC *.cpp)
add_library(ImmersedBoundary SHARED ${ImmersedBoundary_SRC})
add_dependencies(ImmersedBoundary EspressoConfig)

if(CUDA)
  file(GLOB ImmersedBoundaryCuda_SRC *.cu)
  cuda_add_library(ImmersedBoundaryCuda SHARED ${ImmersedBoundaryCuda_SRC})

  set_target_properties(ImmersedBoundaryCuda PROPERTIES MACOSX_RPATH TRUE)
  CUDA_ADD_CUFFT_TO_TARGET(ImmersedBoundaryCuda)
  target_link_libraries(ImmersedBoundary ImmersedBoundaryCuda)
endif(CUDA)

set_target_properties(ImmersedBoundary PROPERTIES LINKER_LANGUAGE CXX)
