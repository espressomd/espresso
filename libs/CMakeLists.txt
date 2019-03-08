add_library(Random123 INTERFACE)
target_include_directories(Random123 SYSTEM INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/Random123-1.09/include)

if(CMAKE_SIZEOF_VOID_P EQUAL 4)
  target_compile_definitions(Random123 INTERFACE R123_USE_MULHILO64_C99)
endif()
