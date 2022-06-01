set(SOURCE_FILES src/Scafacos.cpp src/Coulomb.cpp src/Dipoles.cpp)

add_library(EspressoScafacos SHARED ${SOURCE_FILES})
target_link_libraries(EspressoScafacos PRIVATE ${SCAFACOS_LDFLAGS})
target_link_libraries(EspressoScafacos PUBLIC MPI::MPI_CXX
                      PRIVATE cxx_interface)

target_include_directories(
  EspressoScafacos PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
                          $<INSTALL_INTERFACE:include>)
target_include_directories(EspressoScafacos SYSTEM
                           PUBLIC ${SCAFACOS_INCLUDE_DIRS})

install(TARGETS EspressoScafacos DESTINATION ${PYTHON_INSTDIR}/espressomd)
