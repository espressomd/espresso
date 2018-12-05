add_library(bonded_interactions SHARED
    angle_cosine.cpp
    angle_cossquare.cpp
    angle_harmonic.cpp
    bonded_coulomb.cpp
    bonded_coulomb_p3m_sr.cpp
    bonded_interaction_data.cpp
    bonded_tab.cpp
    dihedral.cpp
    fene.cpp
    harmonic.cpp
    harmonic_dumbbell.cpp
    quartic.cpp
    subt_lj.cpp
    thermalized_bond.cpp
    umbrella.cpp
)

target_link_libraries(bonded_interactions PRIVATE EspressoConfig)
target_include_directories(bonded_interactions PRIVATE ${CMAKE_SOURCE_DIR}/src/core ${CMAKE_BINARY_DIR}/src/core)
