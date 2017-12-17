function(EXAMPLE)
  # Parse arguments
  cmake_parse_arguments(EX "" "TARGET" "SOURCE" ${ARGN} )
  # Add target for test and link with libraries
  add_executable(${EX_TARGET} ${EX_SOURCE})
  if(MSVC)
    set_source_files_properties(${EX_SOURCE} PROPERTIES COMPILE_FLAGS "/W1 /DNOMINMAX")
  else()
    set_source_files_properties(${EX_SOURCE} PROPERTIES COMPILE_FLAGS "-Wall -Wextra -Wpedantic")
  endif()
endfunction(EXAMPLE)

EXAMPLE(TARGET calculator SOURCE calculator.cpp)
EXAMPLE(TARGET motivation SOURCE motivation.cpp)

find_package(Readline)
if(READLINE_FOUND)
  EXAMPLE(TARGET calc SOURCE calc.cpp)
  target_link_libraries(calc ${READLINE_LIBRARY})
endif()
