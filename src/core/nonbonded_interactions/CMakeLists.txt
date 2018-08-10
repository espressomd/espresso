add_library(nonbonded_interactions SHARED
    bmhtf-nacl.cpp
    buckingham.cpp
    cos2.cpp
    debye_hueckel.cpp
    gaussian.cpp
    gb.cpp
    hat.cpp
    hertzian.cpp
    ljcos2.cpp
    ljcos.cpp
    lj.cpp
    ljgen.cpp
    morse.cpp
    nonbonded_interaction_data.cpp
    nonbonded_tab.cpp
    reaction_field.cpp
    soft_sphere.cpp
    steppot.cpp
    thole.cpp
)

target_include_directories(nonbonded_interactions PRIVATE ${CMAKE_SOURCE_DIR}/src/core ${CMAKE_BINARY_DIR}/src/core)
