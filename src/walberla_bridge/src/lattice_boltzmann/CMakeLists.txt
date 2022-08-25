add_subdirectory(generated_kernels)

target_sources(
        espresso_walberla
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/lb_walberla_init.cpp)
