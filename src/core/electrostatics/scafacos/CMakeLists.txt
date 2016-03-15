include_directories(${SCAFACOS_INCLUDE_DIRS})

file(GLOB Scafacos_SRC *.cpp)
add_library(Scafacos SHARED ${Scafacos_SRC})
add_dependencies(Scafacos EspressoConfig)

target_link_libraries(Scafacos ${SCAFACOS_LDFLAGS})
