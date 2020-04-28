set(Scafacos_SRC Scafacos.cpp)
add_library(Scafacos SHARED ${Scafacos_SRC})
target_include_directories(Scafacos PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
                                           ${SCAFACOS_INCLUDE_DIRS})
target_link_libraries(Scafacos PRIVATE ${SCAFACOS_LDFLAGS})
target_link_libraries(Scafacos PUBLIC MPI::MPI_CXX PRIVATE cxx_interface)

install(TARGETS Scafacos DESTINATION ${PYTHON_INSTDIR}/espressomd)
