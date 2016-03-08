add_executable(Vector_test Vector_test.cpp)

add_test(vector_test Vector_test)

add_custom_target(check_unit_tests COMMAND ${CMAKE_CTEST_COMMAND})
add_dependencies(check_unit_tests Vector_test)
add_dependencies(check check_unit_tests)
