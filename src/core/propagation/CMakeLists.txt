target_sources(
  EspressoCore
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/velocity_verlet_npt.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/steepest_descent.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/virtual_sites_relative.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/lb_inertialess_tracers.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/lb_inertialess_tracers_cuda_interface.cpp)
