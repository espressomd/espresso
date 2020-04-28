target_sources(
  EspressoCore
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/lb_inertialess_tracers.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/lb_inertialess_tracers_cuda_interface.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/VirtualSitesInertialessTracers.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/VirtualSitesRelative.cpp)
