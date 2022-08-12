target_sources(
  espresso_script_interface
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/initialize.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/BrownianDynamics.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/IntegratorHandle.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/SteepestDescent.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/StokesianDynamics.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/VelocityVerlet.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/VelocityVerletIsoNPT.cpp)
