target_sources(
  EspressoCore
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/ActorList.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/DipolarBarnesHut.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/DipolarDirectSum.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/Mmm1dgpuForce.cpp)
