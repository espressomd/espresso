file(GLOB ImmersedBoundary_SRC *.cpp)
add_library(ImmersedBoundary SHARED ${ImmersedBoundary_SRC})
set_target_properties(ImmersedBoundary PROPERTIES MACOSX_RPATH TRUE)
install(TARGETS ImmersedBoundary LIBRARY DESTINATION ${PYTHON_INSTDIR} ARCHIVE DESTINATION ${PYTHON_INSTDIR})
add_dependencies(ImmersedBoundary EspressoConfig)

if(CUDA)
  file(GLOB ImmersedBoundaryCuda_SRC *.cu)
  cuda_add_library(ImmersedBoundaryCuda SHARED ${ImmersedBoundaryCuda_SRC})
  install(TARGETS ImmersedBoundaryCuda DESTINATION ${PYTHON_INSTDIR})
  add_dependencies(ImmersedBoundaryCuda EspressoConfig)

  set_target_properties(ImmersedBoundaryCuda PROPERTIES MACOSX_RPATH TRUE)
  CUDA_ADD_CUFFT_TO_TARGET(ImmersedBoundaryCuda)
  target_link_libraries(ImmersedBoundary ImmersedBoundaryCuda)
endif(CUDA)

set_target_properties(ImmersedBoundary PROPERTIES LINKER_LANGUAGE CXX)
