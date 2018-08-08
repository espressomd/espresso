set(Scafacos_SRC
    Scafacos.cpp
)
add_library(Scafacos SHARED ${Scafacos_SRC})
target_include_directories(Scafacos ${SCAFACOS_INCLUDE_DIRS})
install(TARGETS Scafacos DESTINATION ${PYTHON_INSTDIR}/espressomd)
add_dependencies(Scafacos EspressoConfig)

target_link_libraries(Scafacos ${SCAFACOS_LDFLAGS})
