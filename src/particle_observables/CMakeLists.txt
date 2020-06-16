add_library(EspressoParticleObservables INTERFACE)
target_include_directories(
  EspressoParticleObservables SYSTEM
  INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
            $<INSTALL_INTERFACE:include>)

install(TARGETS EspressoParticleObservables
        LIBRARY DESTINATION ${PYTHON_INSTDIR}/espressomd)

if(WITH_TESTS)
  add_subdirectory(tests)
endif()
