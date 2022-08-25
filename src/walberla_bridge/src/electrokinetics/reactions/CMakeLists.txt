add_subdirectory(generated_kernels)

target_sources(
        espresso_walberla
        PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/EKReactionImplBulk.cpp
                ${CMAKE_CURRENT_SOURCE_DIR}/EKReactionImplIndexed.cpp)
set_source_files_properties(reactions/EKReactionImplIndexed.cpp
        TARGET_DIRECTORY espresso_walberla
        PROPERTIES COMPILE_FLAGS -Wno-unused-variable)
