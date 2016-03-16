set(py_tests  bondedInteractions.py  
              electrostaticInteractions.py  
              nonBondedInteractions.py
              cellsystem.py          
              ewald_gpu.py                  
              p3m_gpu.py
              magnetostaticInteractions.py  
              particle.py
)

foreach(testfile ${py_tests})
  get_filename_component(basename ${testfile} NAME_WE)
  add_test(${basename} ${CMAKE_BINARY_DIR}/pypresso ${testfile})
  set_tests_properties(${basename} PROPERTIES SKIP_RETURN_CODE 42)
endforeach(testfile ${py_tests})

add_custom_target(check_python COMMAND ${CMAKE_CTEST_COMMAND})
add_dependencies(check_python pypresso testsuite)
add_dependencies(check check_python)
