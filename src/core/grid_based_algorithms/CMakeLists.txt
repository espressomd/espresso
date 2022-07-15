target_sources(
  Espresso_core
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/lb_interface.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/lb_interpolation.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/lb_particle_coupling.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/lb_walberla_instance.cpp)
