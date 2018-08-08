add_library(H5mdCore SHARED h5md_core.cpp)
target_include_directories(H5mdCore SYSTEM PRIVATE
  # include the h5xx C++ library for h5md
  "${CMAKE_SOURCE_DIR}/libs/h5xx"
  )
add_dependencies(H5mdCore EspressoConfig)

