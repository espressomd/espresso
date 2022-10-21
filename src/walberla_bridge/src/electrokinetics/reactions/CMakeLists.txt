add_subdirectory(generated_kernels)

target_sources(
        espresso_walberla
        PRIVATE EKReactionImplBulk.cpp EKReactionImplIndexed.cpp)
set_source_files_properties(
        EKReactionImplIndexed.cpp
        TARGET_DIRECTORY espresso_walberla
        PROPERTIES COMPILE_FLAGS -Wno-unused-variable)
