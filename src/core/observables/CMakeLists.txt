set(Observables_SRC
    CylindricalLBFluxDensityProfileAtParticlePositions.cpp
    CylindricalLBVelocityProfileAtParticlePositions.cpp
    CylindricalLBVelocityProfile.cpp
    LBVelocityProfile.cpp
)
add_library(Observables SHARED ${Observables_SRC})
install(TARGETS Observables LIBRARY DESTINATION ${PYTHON_INSTDIR}/espressomd)
set_target_properties(Observables PROPERTIES SOVERSION ${SOVERSION})
target_include_directories(Observables PRIVATE ${CMAKE_SOURCE_DIR}/src/core ${CMAKE_BINARY_DIR}/src/core)
target_link_libraries(Observables PRIVATE EspressoConfig)

