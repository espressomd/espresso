add_library(Random123 INTERFACE)
target_include_directories(Random123 SYSTEM INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/Random123-1.09/include)
target_compile_definitions(Random123 INTERFACE R123_USE_MULHILO64_C99)

add_library(h5xx INTERFACE)
target_include_directories(h5xx SYSTEM INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/h5xx)

if(STOKESIAN_DYNAMICS OR STOKESIAN_DYNAMICS_GPU)
  add_library(stokesian_dynamics INTERFACE)
  target_include_directories(stokesian_dynamics INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/stokesian_dynamics/include>)
  set(CMAKE_INSTALL_LIBDIR "${CMAKE_INSTALL_PREFIX}/${PYTHON_INSTDIR}/espressomd")
  add_subdirectory(stokesian_dynamics/src)
endif()
