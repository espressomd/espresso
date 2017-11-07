function(EXAMPLE)
  # Parse arguments
  cmake_parse_arguments(EX "" "TARGET" "SOURCE" ${ARGN} )
  # Add target for test and link with libraries
  add_executable(${EX_TARGET} ${EX_SOURCE})
  set_source_files_properties(${EX_SOURCE} PROPERTIES COMPILE_FLAGS "-Wall -Wextra -Wpedantic")
endfunction(EXAMPLE)

EXAMPLE(TARGET calculator SOURCE calculator.cpp)
EXAMPLE(TARGET motivation SOURCE motivation.cpp)
