foreach(n_species RANGE 1 5)
    foreach(precision double_precision single_precision)
        target_sources(
                espresso_walberla
                PRIVATE
                ReactionKernelBulk_${n_species}_${precision}.cpp
                ReactionKernelIndexed_${n_species}_${precision}.cpp
        )
        set_source_files_properties(
                ReactionKernelBulk_${n_species}_${precision}.cpp
                ReactionKernelIndexed_${n_species}_${precision}.cpp
                PROPERTIES COMPILE_FLAGS -Wno-unused-variable)
    endforeach()
endforeach()