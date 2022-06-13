add_library(Espresso_particle_observables INTERFACE)
add_library(Espresso::particle_observables ALIAS Espresso_particle_observables)
target_include_directories(
  Espresso_particle_observables SYSTEM
  INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
            $<INSTALL_INTERFACE:include>)

install(TARGETS Espresso_particle_observables
        LIBRARY DESTINATION ${PYTHON_INSTDIR}/espressomd)

if(WITH_TESTS)
  add_subdirectory(tests)
endif()
