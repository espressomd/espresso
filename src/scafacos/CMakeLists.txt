add_library(espresso_scafacos SHARED src/Scafacos.cpp src/Coulomb.cpp
                                     src/Dipoles.cpp)
add_library(espresso::scafacos ALIAS espresso_scafacos)
if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 10 AND INSIDE_DOCKER
   AND CMAKE_CXX_COMPILER_VERSION VERSION_LESS_EQUAL 11 AND WITH_COVERAGE)
  target_link_libraries(espresso_scafacos
                        PRIVATE "-L/usr/lib/gcc/x86_64-linux-gnu/10")
endif()
target_link_libraries(espresso_scafacos PUBLIC MPI::MPI_CXX
                      PRIVATE ${SCAFACOS_LDFLAGS} espresso::cpp_flags)

target_include_directories(
  espresso_scafacos
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
         $<INSTALL_INTERFACE:include>)
target_include_directories(espresso_scafacos SYSTEM
                           PUBLIC ${SCAFACOS_INCLUDE_DIRS})

install(TARGETS espresso_scafacos
        DESTINATION ${ESPRESSO_INSTALL_PYTHON}/espressomd)
