include_directories(${SCAFACOS_INCLUDE_DIRS})
add_dependencies(Scafacos EspressoConfig)

file(GLOB Scafacos_SRC *.cpp)
add_library(Scafacos SHARED ${Scafacos_SRC})

target_link_libraries(Scafacos ${SCAFACOS_LDFLAGS})
