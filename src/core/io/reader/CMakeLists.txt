add_library(pdbreader SHARED readpdb.cpp)
target_link_libraries(pdbreader EspressoConfig)
target_include_directories(pdbreader PRIVATE ${CMAKE_SOURCE_DIR}/src/core ${CMAKE_BINARY_DIR}/src/core)
