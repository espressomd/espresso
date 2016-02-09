file(GLOB all_files *)
add_custom_target(testsuite_py ALL)
foreach(testfile ${all_files})
  get_filename_component(basename ${testfile} NAME)
  add_custom_command(TARGET testsuite_py
                     COMMAND test -d ${basename} || ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_SOURCE_DIR}/${basename} ${basename})
endforeach(testfile)

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
  add_test(${basename} ../../pypresso ${testfile})
  set_tests_properties(${basename} PROPERTIES SKIP_RETURN_CODE 42)
endforeach(testfile ${py_tests})
