add_library(pdbparser SHARED PdbParser.cpp)
target_include_directories(pdbparser PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
