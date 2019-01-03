set(ImmersedBoundary_SRC
    ibm_tribend.cpp
    ibm_triel.cpp
    ImmersedBoundaries.cpp
)
add_library(ImmersedBoundary SHARED ${ImmersedBoundary_SRC})
install(TARGETS ImmersedBoundary LIBRARY DESTINATION ${PYTHON_INSTDIR}/espressomd ARCHIVE DESTINATION ${PYTHON_INSTDIR}/espressomd)
set_target_properties(ImmersedBoundary PROPERTIES SOVERSION ${SOVERSION})
add_dependencies(ImmersedBoundary EspressoConfig)
set_target_properties(ImmersedBoundary PROPERTIES LINKER_LANGUAGE CXX)
target_include_directories(ImmersedBoundary PRIVATE ${CMAKE_SOURCE_DIR}/src/core ${CMAKE_BINARY_DIR}/src/core)
