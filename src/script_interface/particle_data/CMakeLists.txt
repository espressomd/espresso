target_sources(
  Espresso_script_interface
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/initialize.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/ParticleHandle.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/ParticleList.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/ParticleSlice.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/Polymer.cpp)
