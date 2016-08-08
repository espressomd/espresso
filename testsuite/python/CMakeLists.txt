set(py_tests  bondedInteractions.py
              cellsystem.py
              coulomb_cloud_wall.py
              electrostaticInteractions.py
              engine_langevin.py
              engine_lb.py
              ewald_gpu.py
              magnetostaticInteractions.py
              nonBondedInteractions.py
              p3m_gpu.py
              particle.py
)

add_custom_target(python_tests
                  COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})

foreach(testfile ${py_tests})
  get_filename_component(basename ${testfile} NAME_WE)
  add_test(${basename} ${CMAKE_BINARY_DIR}/pypresso ${testfile})
endforeach(testfile ${py_tests})

add_custom_target(check_python COMMAND ${CMAKE_CTEST_COMMAND})
add_dependencies(check_python pypresso python_tests)
add_dependencies(check check_python)
