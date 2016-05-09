add_executable(Vector_test Vector_test.cpp)
add_executable(RuntimeError_test RuntimeError_test.cpp)
add_executable(RuntimeErrorCollector_test RuntimeErrorCollector_test.cpp ../RuntimeErrorCollector.cpp ../RuntimeError.cpp)
target_link_libraries(RuntimeErrorCollector_test ${LIBRARIES})
target_link_libraries(RuntimeError_test ${LIBRARIES})

add_test(vector_test Vector_test)
add_test(runtimeerror_test RuntimeError_test)
add_test(runtimeerrorcollector_test RuntimeErrorCollector_test)

add_custom_target(check_unit_tests COMMAND ${CMAKE_CTEST_COMMAND})
add_dependencies(check_unit_tests Vector_test)
add_dependencies(check_unit_tests RuntimeError_test)
add_dependencies(check_unit_tests RuntimeErrorCollector_test)

add_dependencies(check check_unit_tests)
