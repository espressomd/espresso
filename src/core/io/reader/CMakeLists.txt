add_library(pdbreader SHARED readpdb.cpp)
target_link_libraries(pdbreader PUBLIC EspressoConfig EspressoCore pdbparser)
target_include_directories(pdbreader PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
