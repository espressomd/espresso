if (HDF5_FOUND)
add_library(H5mdCore SHARED h5md_core.cpp)
target_include_directories(H5mdCore SYSTEM PUBLIC
  # include the h5xx C++ library for h5md
  "${CMAKE_SOURCE_DIR}/libs/h5xx"
)

target_include_directories(H5mdCore PUBLIC SYSTEM ${HDF5_INCLUDE_DIRS})
target_link_libraries(H5mdCore PUBLIC EspressoConfig EspressoCore ${HDF5_LIBRARIES} Boost::filesystem)
install(TARGETS H5mdCore LIBRARY DESTINATION ${PYTHON_INSTDIR}/espressomd)
endif()
