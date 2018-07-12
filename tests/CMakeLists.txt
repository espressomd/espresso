# Build tests only when testing
set_directory_properties(PROPERTIES EXCLUDE_FROM_ALL ON)

# Add current directory to search path
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# To add test cases easier we use the following function
function(UNIT_TEST)
  # Parse arguments
  cmake_parse_arguments(UT "" "TARGET" "SOURCE" ${ARGN} )
  # Add target for test and link with libraries
  add_executable(${UT_TARGET} ${UT_SOURCE})
  set_target_properties(${UT_TARGET} PROPERTIES CXX_CLANG_TIDY "")
  if(MSVC)
    set_property(SOURCE ${UT_SOURCE} PROPERTY COMPILE_FLAGS "/DNOMINMAX")
  endif()
  target_link_libraries(${UT_TARGET} ${LIBRARIES} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
  # Tests are executed in the root directory
  add_test(NAME ${UT_TARGET} COMMAND ${UT_TARGET})
  add_dependencies(check ${UT_TARGET})
endfunction(UNIT_TEST)

# Add the tests
unit_test(TARGET literals SOURCE literals.cpp)
unit_test(TARGET arithmetic SOURCE arithmetic.cpp)
unit_test(TARGET order_of_operations SOURCE order_of_operations.cpp)
unit_test(TARGET unary SOURCE unary.cpp)
unit_test(TARGET functions SOURCE functions.cpp)
unit_test(TARGET variables SOURCE variables.cpp)
unit_test(TARGET optimizer SOURCE optimizer.cpp)
unit_test(TARGET errors SOURCE errors.cpp)
unit_test(TARGET empty_node SOURCE empty_node.cpp)
unit_test(TARGET interface SOURCE interface.cpp)
