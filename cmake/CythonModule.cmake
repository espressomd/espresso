function(ADD_CYTHON_MODULE)
  cmake_parse_arguments(CY "" "PYX;TARGET_DIR;TARGET_PREFIX" "DEPENDS;LINK_WITH" ${ARGN})
  get_filename_component(basename ${CY_PYX} NAME_WE)
  set(outputpath "${CY_TARGET_DIR}/${basename}.cpp")
  add_custom_command(OUTPUT ${outputpath}
                     COMMAND ${CYTHON_EXECUTABLE} --cplus
                     --directive embedsignature=True
                     --directive binding=True
                     -I ${CMAKE_CURRENT_SOURCE_DIR}
                     -I ${CMAKE_CURRENT_BINARY_DIR}
                     ${CY_PYX}
                     -o ${outputpath}
                     DEPENDS ${CY_DEPENDS} ${CY_PYX}
                     )
  set(target "${CY_TARGET_PREFIX}_${basename}")
  add_library(${target} SHARED ${outputpath})
  set_target_properties(${target} PROPERTIES OUTPUT_NAME ${basename})
  if (APPLE)
     set_target_properties(${target} PROPERTIES
             SUFFIX ".so"
             LINK_FLAGS "-undefined dynamic_lookup"
             )
  endif()
  if ( CMAKE_CXX_COMPILER_ID STREQUAL "GNU" )
      set_source_files_properties(${outputpath} PROPERTIES COMPILE_FLAGS "-Wno-pedantic -Wno-cpp -Wno-strict-aliasing -Wno-maybe-uninitialized -Wno-unused-variable")
  elseif( CMAKE_CXX_COMPILER_ID STREQUAL "Clang" OR CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang" )
      set_source_files_properties(${outputpath} PROPERTIES COMPILE_FLAGS "-Wno-pedantic -Wno-#warnings -Wno-sometimes-uninitialized -Wno-unused-variable")
  elseif( CMAKE_CXX_COMPILER_ID STREQUAL "Intel" )
      set_source_files_properties(${outputpath} PROPERTIES COMPILE_FLAGS "-wd1224")
  else()
      set_source_files_properties(${outputpath} PROPERTIES COMPILE_FLAGS "-Wno-pedantic -Wno-unused-variable")
  endif()
  set_target_properties(${target} PROPERTIES CXX_CLANG_TIDY "")
  target_link_libraries(${target} PRIVATE ${CY_LINK_WITH})
  target_include_directories(${target} SYSTEM PRIVATE ${PYTHON_INCLUDE_DIRS} ${NUMPY_INCLUDE_DIR})
  add_dependencies(espressomd ${target})
endfunction()
