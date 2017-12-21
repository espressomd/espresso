file(GLOB ImmersedBoundary_SRC *.cpp)
add_library(ImmersedBoundary SHARED ${ImmersedBoundary_SRC})
install(TARGETS ImmersedBoundary LIBRARY DESTINATION ${LIBDIR} ARCHIVE DESTINATION ${LIBDIR})
set_target_properties(ImmersedBoundary PROPERTIES SOVERSION ${SOVERSION})
add_dependencies(ImmersedBoundary EspressoConfig)

if(CUDA)
  file(GLOB ImmersedBoundaryCuda_SRC *.cu)
  cuda_add_library(ImmersedBoundaryCuda SHARED ${ImmersedBoundaryCuda_SRC})
  install(TARGETS ImmersedBoundaryCuda DESTINATION ${LIBDIR})
  set_target_properties(ImmersedBoundaryCuda PROPERTIES SOVERSION ${SOVERSION})
  add_dependencies(ImmersedBoundaryCuda EspressoConfig)

  CUDA_ADD_CUFFT_TO_TARGET(ImmersedBoundaryCuda)
  target_link_libraries(ImmersedBoundary ImmersedBoundaryCuda)
endif(CUDA)

set_target_properties(ImmersedBoundary PROPERTIES LINKER_LANGUAGE CXX)
