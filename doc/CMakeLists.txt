add_subdirectory(doxygen)
add_subdirectory(sphinx)
add_subdirectory(logo)
add_subdirectory(tutorials)

add_custom_target(doc DEPENDS logo doxygen tutorials)
