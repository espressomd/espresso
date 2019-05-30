function(TUTORIAL_TEST)
  PYTHON_SCRIPTS_TEST(${ARGV} TYPE tutorial)  # forward arguments
  set(tutorials_tests ${tutorials_tests} ${TEST_FILE_CONFIGURED} PARENT_SCOPE)
endfunction(TUTORIAL_TEST)

# configure Python module importlib_wrapper.py
add_test(importlib_wrapper ${CMAKE_BINARY_DIR}/pypresso ${TEST_FILE_CONFIGURED_IMPORTLIB_WRAPPER})
set_tests_properties(importlib_wrapper PROPERTIES FIXTURES_SETUP IMPORTLIB_WRAPPER)
set(tutorials_tests ${tutorials_tests} ${TEST_FILE_CONFIGURED_IMPORTLIB_WRAPPER} PARENT_SCOPE)
configure_file(../importlib_wrapper.py ${CMAKE_CURRENT_BINARY_DIR}/importlib_wrapper.py)

# create a local copy of the tutorials
set(TUTORIALS_DIR "${CMAKE_CURRENT_BINARY_DIR}/local_tutorials")
add_custom_target(local_tutorials
                  COMMAND ${CMAKE_COMMAND} -E remove_directory ${TUTORIALS_DIR}  # cleanup
                  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_BINARY_DIR}/doc/tutorials ${TUTORIALS_DIR}
                  DEPENDS tutorials_python_only)

tutorial_test(FILE test_01-lennard_jones.py)
tutorial_test(FILE test_02-charged_system-1.py)
tutorial_test(FILE test_02-charged_system-2.py)
tutorial_test(FILE test_02-charged_system__scripts__nacl.py)
tutorial_test(FILE test_02-charged_system__scripts__nacl_units_confined.py)
tutorial_test(FILE test_02-charged_system__scripts__nacl_units_confined_vis.py)
tutorial_test(FILE test_02-charged_system__scripts__nacl_units.py)
tutorial_test(FILE test_02-charged_system__scripts__nacl_units_vis.py)
tutorial_test(FILE test_04-lattice_boltzmann_part1.py LABELS "gpu")
tutorial_test(FILE test_04-lattice_boltzmann_part2.py LABELS "gpu")
tutorial_test(FILE test_04-lattice_boltzmann_part3.py LABELS "gpu")
tutorial_test(FILE test_04-lattice_boltzmann_part4.py LABELS "gpu")
tutorial_test(FILE test_05-raspberry_electrophoresis.py LABELS "gpu")
tutorial_test(FILE test_06-active_matter__flow_field.py LABELS "gpu")
tutorial_test(FILE test_06-active_matter__rectification_geometry.py LABELS "gpu")
tutorial_test(FILE test_06-active_matter__enhanced_diffusion.py)
tutorial_test(FILE test_06-active_matter__rectification_simulation.py)
tutorial_test(FILE test_07-electrokinetics__scripts.py LABELS "gpu")
tutorial_test(FILE test_08-visualization.py)
tutorial_test(FILE test_09-swimmer_reactions__reaction.py)
tutorial_test(FILE test_10-reaction_ensemble__scripts__RE_vs_cpH.py)
tutorial_test(FILE test_10-reaction_ensemble__scripts__RE_vs_cpH_poly.py)

add_custom_target(check_tutorials COMMAND ${CMAKE_CTEST_COMMAND} --timeout ${TEST_TIMEOUT} -j${TEST_NP} $(ARGS) --output-on-failure)

add_dependencies(check_tutorials pypresso local_tutorials)

