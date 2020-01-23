add_custom_target(examples)
add_dependencies(examples matheval.qi matheval.x3)

function(EXAMPLE)
  # Parse arguments
  cmake_parse_arguments(EX "" "TARGET" "SOURCE" ${ARGN} )
  # Add target for test and link with libraries
  add_executable(matheval.ex.${EX_TARGET} ${EX_SOURCE})
  target_link_libraries(matheval.ex.${EX_TARGET} PRIVATE matheval::x3)
  set_target_properties(matheval.ex.${EX_TARGET} PROPERTIES CXX_STANDARD 14)
  if(MSVC)
    set_source_files_properties(${EX_SOURCE} PROPERTIES COMPILE_FLAGS "/W1 /DNOMINMAX")
  else()
    set_source_files_properties(${EX_SOURCE} PROPERTIES COMPILE_FLAGS "-Wall -Wextra -Wpedantic")
  endif()
  add_dependencies(examples matheval.ex.${EX_TARGET})
endfunction(EXAMPLE)

EXAMPLE(TARGET calculator SOURCE calculator.cpp)
EXAMPLE(TARGET motivation SOURCE motivation.cpp)

find_package(Readline)
if(READLINE_FOUND)
  EXAMPLE(TARGET calc SOURCE calc.cpp)
  target_link_libraries(matheval.ex.calc PRIVATE ${READLINE_LIBRARY})
endif()
