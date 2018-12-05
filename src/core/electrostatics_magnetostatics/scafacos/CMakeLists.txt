set(Scafacos_SRC
    Scafacos.cpp
)
add_library(Scafacos SHARED ${Scafacos_SRC})
if (${SCAFACOS_INCLUDE_DIRS})
    target_include_directories(Scafacos PUBLIC ${SCAFACOS_INCLUDE_DIRS})
endif()
install(TARGETS Scafacos DESTINATION ${PYTHON_INSTDIR}/espressomd)
target_link_libraries(Scafacos PRIVATE ${SCAFACOS_LDFLAGS})
