configure_file("pypresso.cmakein" "pypresso")
add_custom_target(pypresso ALL
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/pypresso)
add_custom_command(TARGET pypresso
  COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_BINARY_DIR}/pypresso ${CMAKE_BINARY_DIR}/pypresso
  COMMAND ${CMAKE_COMMAND} -E create_symlink ${CMAKE_CURRENT_BINARY_DIR}/espressomd ${CMAKE_BINARY_DIR}/espressomd)
add_subdirectory(espressomd)
