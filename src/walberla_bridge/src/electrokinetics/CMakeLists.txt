add_subdirectory(generated_kernels)
add_subdirectory(reactions)

target_sources(
        espresso_walberla
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/ek_walberla_init.cpp)

set_source_files_properties(
        src/electrokinetics/ek_walberla_init.cpp
        PROPERTIES COMPILE_FLAGS -Wno-unused-variable)