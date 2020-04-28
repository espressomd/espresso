target_sources(
  EspressoCore
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/halo.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/lattice.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/lb_boundaries.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/lb_collective_interface.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/lb.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/lb_interface.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/lb_interpolation.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/lb_particle_coupling.cpp)
