target_include_directories(
  EspressoCore PRIVATE ${CMAKE_SOURCE_DIR}/libs
                       ${CMAKE_SOURCE_DIR}/src/core/stokesian_dynamics)

target_sources(EspressoCore
               PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/sd_interface.cpp)
