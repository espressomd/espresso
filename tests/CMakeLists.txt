# Build tests only when testing
set_directory_properties(PROPERTIES EXCLUDE_FROM_ALL ON)

# Add current directory to search path
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# To add test cases easier we use the following function
function(UNIT_TEST)
  # Parse arguments
  cmake_parse_arguments(UT "" "TARGET" "SOURCE" ${ARGN} )
  # Add target for QI
  add_executable(qi.${UT_TARGET} ${UT_SOURCE})
  target_include_directories(qi.${UT_TARGET} PRIVATE ${CMAKE_SOURCE_DIR}/include/matheval/qi)
  target_link_libraries(qi.${UT_TARGET} PRIVATE matheval_qi)
  set_target_properties(qi.${UT_TARGET} PROPERTIES CXX_CLANG_TIDY "" CXX_STANDARD 98)
  if(MSVC)
    set_property(SOURCE ${UT_SOURCE} PROPERTY COMPILE_FLAGS "/DNOMINMAX")
  endif()
  # Add target for X3
  add_executable(x3.${UT_TARGET} ${UT_SOURCE})
  target_include_directories(x3.${UT_TARGET} PRIVATE ${CMAKE_SOURCE_DIR}/include/matheval/x3)
  target_link_libraries(x3.${UT_TARGET} PRIVATE matheval_x3)
  set_target_properties(x3.${UT_TARGET} PROPERTIES CXX_CLANG_TIDY "" CXX_STANDARD 14)
  if(MSVC)
    set_property(SOURCE ${UT_SOURCE} PROPERTY COMPILE_FLAGS "/DNOMINMAX")
  endif()
  # Tests are executed in the root directory
  add_test(NAME qi.${UT_TARGET} COMMAND qi.${UT_TARGET})
  add_test(NAME x3.${UT_TARGET} COMMAND x3.${UT_TARGET})
  add_dependencies(check qi.${UT_TARGET} x3.${UT_TARGET})
endfunction(UNIT_TEST)

# Add the tests
unit_test(TARGET literals SOURCE literals.cpp)
unit_test(TARGET arithmetic SOURCE arithmetic.cpp)
unit_test(TARGET order_of_operations SOURCE order_of_operations.cpp)
unit_test(TARGET unary SOURCE unary.cpp)
unit_test(TARGET functions SOURCE functions.cpp)
unit_test(TARGET variables SOURCE variables.cpp)
#unit_test(TARGET optimizer SOURCE optimizer.cpp)
unit_test(TARGET errors SOURCE errors.cpp)
#unit_test(TARGET empty_node SOURCE empty_node.cpp)
unit_test(TARGET interface SOURCE interface.cpp)
