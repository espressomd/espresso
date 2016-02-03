file(GLOB Actor_SRC *.cpp)

add_library(Actor SHARED ${Actor_SRC})

if(WITH_CUDA)
  file(GLOB ActorCuda_SRC *.cu)
  cuda_add_library(ActorCuda SHARED ${ActorCuda_SRC})

  add_dependencies(Actor ActorCuda)
  set_target_properties(ActorCuda PROPERTIES MACOSX_RPATH TRUE)
  CUDA_ADD_CUFFT_TO_TARGET(ActorCuda)
  target_link_libraries(Actor ActorCuda)
endif(WITH_CUDA)
