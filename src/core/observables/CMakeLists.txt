set(Observables_SRC
    CylindricalLBFluxDensityProfileAtParticlePositions.cpp
    CylindricalLBVelocityProfileAtParticlePositions.cpp
    CylindricalLBVelocityProfile.cpp
    LBVelocityProfile.cpp
)
add_library(Observables SHARED ${Observables_SRC})
install(TARGETS Observables LIBRARY DESTINATION ${LIBDIR})
set_target_properties(Observables PROPERTIES SOVERSION ${SOVERSION})
add_dependencies(Observables EspressoConfig)
