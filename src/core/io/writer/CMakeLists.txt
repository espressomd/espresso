if (HDF5_FOUND)
add_library(H5mdCore SHARED h5md_core.cpp)
target_include_directories(H5mdCore SYSTEM PUBLIC
  # include the h5xx C++ library for h5md
  "${CMAKE_SOURCE_DIR}/libs/h5xx"
)

target_include_directories(H5mdCore PRIVATE ${CMAKE_SOURCE_DIR}/src/core ${CMAKE_BINARY_DIR}/src/core)
target_link_libraries(H5mdCore PRIVATE EspressoConfig)
install(TARGETS H5mdCore LIBRARY DESTINATION ${PYTHON_INSTDIR}/espressomd)
endif()
