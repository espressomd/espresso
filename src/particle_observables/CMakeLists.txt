add_library(espresso_particle_observables INTERFACE)
add_library(espresso::particle_observables ALIAS espresso_particle_observables)
target_include_directories(
  espresso_particle_observables SYSTEM
  INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
            $<INSTALL_INTERFACE:include>)

install(TARGETS espresso_particle_observables
        LIBRARY DESTINATION ${ESPRESSO_INSTALL_PYTHON}/espressomd)

if(WITH_TESTS)
  add_subdirectory(tests)
endif()
