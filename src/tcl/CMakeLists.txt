file(GLOB_RECURSE Espresso_SRC
          "*.cpp"
          )

include_directories(${CMAKE_SOURCE_DIR}/src/core)
include_directories(${CMAKE_BINARY_DIR}/src/core)
include_directories(${CMAKE_SOURCE_DIR}/src/tcl)

file(GLOB espresso_scripts ${CMAKE_SOURCE_DIR}/scripts/*tcl)

file(COPY ${espresso_scripts} DESTINATION ${CMAKE_BINARY_DIR}/scripts)
install(DIRECTORY ${CMAKE_BINARY_DIR}/scripts DESTINATION ${DATA})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
add_definitions(-DESPRESSO_SCRIPTS_DEFAULT="${CMAKE_BINARY_DIR}/scripts")
add_executable(Espresso ${Espresso_SRC} ${MPI_SOURCE})
install(TARGETS Espresso RUNTIME DESTINATION bin)
target_link_libraries(Espresso EspressoCore ${LIBRARIES})

