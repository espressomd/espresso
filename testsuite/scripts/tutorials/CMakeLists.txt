function(TUTORIAL_TEST)
  python_scripts_test(${ARGV} TYPE tutorial) # forward arguments
  set(tutorials_tests ${tutorials_tests} ${TEST_FILE_CONFIGURED} PARENT_SCOPE)
endfunction(TUTORIAL_TEST)

# configure Python modules importlib_wrapper.py and html_runner.py
configure_file(test_html_runner.py ${CMAKE_CURRENT_BINARY_DIR})
add_test(importlib_wrapper ${CMAKE_BINARY_DIR}/pypresso
         ${TEST_FILE_CONFIGURED_IMPORTLIB_WRAPPER})
add_test(html_runner ${CMAKE_BINARY_DIR}/pypresso
         ${CMAKE_CURRENT_BINARY_DIR}/test_html_runner.py)
set_tests_properties(importlib_wrapper PROPERTIES FIXTURES_SETUP
                                                  IMPORTLIB_WRAPPER)
set_tests_properties(html_runner PROPERTIES FIXTURES_SETUP IMPORTLIB_WRAPPER)
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

tutorial_test(FILE test_01-lennard_jones.py)
tutorial_test(FILE test_02-charged_system-1.py)
tutorial_test(FILE test_02-charged_system-2.py)
tutorial_test(FILE test_02-charged_system__scripts__nacl.py)
tutorial_test(FILE test_02-charged_system__scripts__nacl_units_confined.py)
tutorial_test(FILE test_02-charged_system__scripts__nacl_units_confined_vis.py)
tutorial_test(FILE test_02-charged_system__scripts__nacl_units.py)
tutorial_test(FILE test_02-charged_system__scripts__nacl_units_vis.py)
tutorial_test(FILE test_04-lattice_boltzmann_part2.py LABELS "gpu")
tutorial_test(FILE test_04-lattice_boltzmann_part3.py LABELS "gpu")
tutorial_test(FILE test_04-lattice_boltzmann_part3_solution.py LABELS "gpu")
tutorial_test(FILE test_04-lattice_boltzmann_part4.py LABELS "gpu")
tutorial_test(FILE test_04-lattice_boltzmann_part4_solution.py LABELS "gpu")
tutorial_test(FILE test_05-raspberry_electrophoresis.py LABELS "gpu")
tutorial_test(FILE test_06-active_matter__flow_field.py LABELS "gpu")
tutorial_test(FILE test_06-active_matter__rectification_geometry.py LABELS
              "gpu")
tutorial_test(FILE test_06-active_matter__enhanced_diffusion.py)
tutorial_test(FILE test_06-active_matter__rectification_simulation.py)
tutorial_test(FILE test_07-electrokinetics.py LABELS "gpu")
tutorial_test(FILE test_08-visualization.py)
tutorial_test(FILE test_11-ferrofluid_1.py)
tutorial_test(FILE test_11-ferrofluid_2.py)
tutorial_test(FILE test_11-ferrofluid_3.py)
tutorial_test(FILE test_12-constant_pH.py)

add_custom_target(
  check_tutorials COMMAND ${CMAKE_CTEST_COMMAND} --timeout ${TEST_TIMEOUT}
                          -j${TEST_NP} $(ARGS) --output-on-failure)

add_dependencies(check_tutorials pypresso local_tutorials)
