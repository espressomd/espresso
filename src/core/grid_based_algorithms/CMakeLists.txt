target_sources(
        EspressoCore
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/lb_boundaries.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/lbboundaries/LBBoundary.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/lb_interface.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/lb_interpolation.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/lb_particle_coupling.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/lb_walberla_interface.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/lb_walberla_instance.cpp
        ${CMAKE_CURRENT_SOURCE_DIR}/ek_container.cpp)
