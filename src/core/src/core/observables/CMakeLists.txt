target_sources(
  EspressoCore
  PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/CylindricalLBFluxDensityProfileAtParticlePositions.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/CylindricalLBVelocityProfileAtParticlePositions.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/CylindricalLBVelocityProfile.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/LBVelocityProfile.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PidObservable.cpp)
