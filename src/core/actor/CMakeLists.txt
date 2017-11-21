file(GLOB Actor_SRC *.cpp)

add_library(Actor SHARED ${Actor_SRC})
install(TARGETS Actor LIBRARY DESTINATION ${LIBDIR} ARCHIVE DESTINATION ${LIBDIR})
set_target_properties(Actor PROPERTIES SOVERSION ${SOVERSION})
add_dependencies(Actor EspressoConfig)

if(CUDA)
  file(GLOB ActorCuda_SRC *.cu)
  cuda_add_library(ActorCuda SHARED ${ActorCuda_SRC})
  install(TARGETS ActorCuda DESTINATION ${LIBDIR})
  set_target_properties(ActorCuda PROPERTIES SOVERSION ${SOVERSION})
  add_dependencies(ActorCuda EspressoConfig)

  add_dependencies(Actor ActorCuda)
  CUDA_ADD_CUFFT_TO_TARGET(ActorCuda)
  target_link_libraries(Actor ActorCuda)
endif(CUDA)
