add_library(EspressoInteractions INTERFACE)
target_include_directories(
  EspressoInteractions
  INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
            $<INSTALL_INTERFACE:include>)

install(TARGETS EspressoInteractions
        LIBRARY DESTINATION ${PYTHON_INSTDIR}/espressomd)

if(WITH_TESTS)
  add_subdirectory(tests)
endif()
