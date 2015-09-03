file(GLOB Espresso_SRC
          "*.cpp"
          )

include_directories(${CMAKE_SOURCE_DIR}/src/core)
include_directories(${CMAKE_BINARY_DIR}/src/core)
include_directories(${CMAKE_SOURCE_DIR}/src/tcl)

add_executable(Espresso ${Espresso_SRC})
target_link_libraries(Espresso libEspresso tcl)
