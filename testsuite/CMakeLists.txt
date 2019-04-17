add_subdirectory(python)
add_subdirectory(cmake)

if(${CMAKE_VERSION} VERSION_LESS "3.7.0")
  add_custom_target(check_tutorials COMMAND echo "Target 'check_tutorials' requires CMake version 3.7.0 or above, currently using ${CMAKE_VERSION}")
  add_custom_target(check_samples COMMAND echo "Target 'check_samples' requires CMake version 3.7.0 or above, currently using ${CMAKE_VERSION}")
else()
  add_subdirectory(scripts)
endif()
