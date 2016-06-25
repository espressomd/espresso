set_property(GLOBAL APPEND PROPERTY LIBRARIES H5mdCore)
add_library(H5mdCore SHARED h5md_core.cpp)
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
