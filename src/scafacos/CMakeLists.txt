set(SOURCE_FILES src/Scafacos.cpp src/Coulomb.cpp src/Dipoles.cpp)

add_library(Espresso_scafacos SHARED ${SOURCE_FILES})
add_library(Espresso::scafacos ALIAS Espresso_scafacos)
target_link_libraries(Espresso_scafacos PRIVATE ${SCAFACOS_LDFLAGS})
target_link_libraries(Espresso_scafacos PUBLIC MPI::MPI_CXX
                      PRIVATE Espresso::cpp_flags)

target_include_directories(
  Espresso_scafacos
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
         $<INSTALL_INTERFACE:include>)
target_include_directories(Espresso_scafacos SYSTEM
                           PUBLIC ${SCAFACOS_INCLUDE_DIRS})

install(TARGETS Espresso_scafacos DESTINATION ${PYTHON_INSTDIR}/espressomd)
