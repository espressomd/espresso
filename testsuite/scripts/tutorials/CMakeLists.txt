function(TUTORIAL_TEST)
  python_scripts_test(${ARGV} TYPE tutorial) # forward arguments
  set(tutorials_tests ${tutorials_tests} ${TEST_FILE_CONFIGURED} PARENT_SCOPE)
endfunction(TUTORIAL_TEST)

# configure Python modules importlib_wrapper.py and convert.py
configure_file(test_convert.py ${CMAKE_CURRENT_BINARY_DIR})
add_test(importlib_wrapper ${CMAKE_BINARY_DIR}/pypresso ${PYPRESSO_OPTIONS}
         ${TEST_FILE_CONFIGURED_IMPORTLIB_WRAPPER})
add_test(convert ${CMAKE_BINARY_DIR}/pypresso ${PYPRESSO_OPTIONS}
         ${CMAKE_CURRENT_BINARY_DIR}/test_convert.py)
set_tests_properties(importlib_wrapper PROPERTIES FIXTURES_SETUP
                                                  IMPORTLIB_WRAPPER)
set_tests_properties(convert PROPERTIES FIXTURES_SETUP IMPORTLIB_WRAPPER)
set(tutorials_tests ${tutorials_tests} ${TEST_FILE_CONFIGURED_IMPORTLIB_WRAPPER}
    PARENT_SCOPE)
configure_file(../importlib_wrapper.py
               ${CMAKE_CURRENT_BINARY_DIR}/importlib_wrapper.py)

# create a local copy of the tutorials
set(TUTORIALS_DIR "${CMAKE_CURRENT_BINARY_DIR}/local_tutorials")
add_custom_target(
  local_tutorials COMMAND ${CMAKE_COMMAND} -E remove_directory
                          ${TUTORIALS_DIR} # cleanup
  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_BINARY_DIR}/doc/tutorials
          ${TUTORIALS_DIR} DEPENDS tutorials_python)

tutorial_test(FILE test_error_analysis_part1.py)
tutorial_test(FILE test_error_analysis_part2.py)
tutorial_test(FILE test_lennard_jones.py)
tutorial_test(FILE test_charged_system.py)
tutorial_test(FILE test_lattice_boltzmann_part2.py)
tutorial_test(FILE test_lattice_boltzmann_part3.py SUFFIX rouse)
tutorial_test(FILE test_lattice_boltzmann_part3.py SUFFIX zimm LABELS "gpu")
tutorial_test(FILE test_lattice_boltzmann_part4.py LABELS "gpu")
tutorial_test(FILE test_raspberry_electrophoresis.py LABELS "gpu")
tutorial_test(FILE test_active_matter__flow_field.py LABELS "gpu")
tutorial_test(FILE test_active_matter__rectification_geometry.py)
tutorial_test(FILE test_active_matter__enhanced_diffusion.py)
tutorial_test(FILE test_active_matter__rectification_simulation.py)
tutorial_test(FILE test_electrokinetics.py LABELS "gpu")
tutorial_test(FILE test_visualization.py)
tutorial_test(FILE test_ferrofluid_1.py)
tutorial_test(FILE test_ferrofluid_2.py)
tutorial_test(FILE test_ferrofluid_3.py)
tutorial_test(FILE test_constant_pH__ideal.py)
tutorial_test(FILE test_constant_pH__interactions.py)

add_custom_target(
  check_tutorials COMMAND ${CMAKE_CTEST_COMMAND} --timeout ${TEST_TIMEOUT}
                          ${CTEST_ARGS} --output-on-failure)

add_dependencies(check_tutorials pypresso local_tutorials)
