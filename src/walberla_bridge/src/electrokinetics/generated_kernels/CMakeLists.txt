foreach(precision double_precision single_precision)
    target_sources(
            espresso_walberla
            PRIVATE
            DiffusiveFluxKernel_${precision}.cpp
            DiffusiveFluxKernelWithElectrostatic_${precision}.cpp
            ContinuityKernel_${precision}.cpp
            AdvectiveFluxKernel_${precision}.cpp
            FrictionCouplingKernel_${precision}.cpp
            FixedFlux_${precision}.cpp
            Dirichlet_${precision}.cpp)
    set_source_files_properties(
            Dirichlet_${precision}.cpp
            Dirichlet_${precision}.h
            PROPERTIES COMPILE_FLAGS -Wno-unused-variable)

endforeach()