enable_testing()
set(CMAKE_TESTS_TO_BE_RUN 
    analysis.tcl
    comforce.tcl
    comfixed.tcl
    command_syntax.tcl
    constraints.tcl
    constraints_reflecting.tcl
    constraints_rhomboid.tcl
    dh.tcl
    el2d.tcl
    el2d_die.tcl
    el2d_nonneutral.tcl
    fene.tcl
    gb.tcl
    harm.tcl
    intpbc.tcl
    intppbc.tcl
    kinetic.tcl
    layered.tcl
    lb.tcl
    lb_stokes_sphere.tcl
    lb_gpu.tcl
    lj.tcl
    lj-cos.tcl
    lj-generic.tcl
    madelung.tcl
    maggs.tcl
    magnetic-field.tcl
    mass.tcl
    mass-and-rinertia.tcl
    mdlc.tcl
    mmm1d.tcl
    npt.tcl
    nsquare.tcl
    nve_pe.tcl
    p3m.tcl
    p3m_magnetostatics.tcl
    p3m_magnetostatics2.tcl
    p3m_simple_noncubic.tcl
    p3m_wall.tcl
    rotation.tcl
    tabulated.tcl
    thermostat.tcl
    tunable_slip.tcl
    virtual-sites.tcl
)

# copy all test files to bin dir
add_custom_target(testfiles COMMENT "Copying Tcl test scripts to builddir")
file(GLOB TESTING_FILES *.gz *.chk *.data *.tab *.tcl)
foreach(FILENAME ${TESTING_FILES})
    add_custom_command(TARGET testfiles 
        COMMAND ${CMAKE_COMMAND} -E copy_if_different ${FILENAME} ${CMAKE_CURRENT_BINARY_DIR}
    )
    install(FILES ${FILENAME} DESTINATION testsuite)
endforeach(FILENAME)

# create runtest.sh file in bin dir
configure_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/runtest.sh.in 
    ${CMAKE_CURRENT_BINARY_DIR}/runtest.sh
)

#add tests
add_custom_target(check
    DEPENDS testfiles Espresso
    COMMAND /bin/sh runtest.sh -p ${MPIEXEC_MAX_NUMPROCS} ${CMAKE_TESTS_TO_BE_RUN}
    COMMENT "Run some ESPResSo test cases on available processors."
)
