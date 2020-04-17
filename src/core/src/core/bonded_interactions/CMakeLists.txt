target_sources(
  EspressoCore
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/angle_cosine.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/angle_cossquare.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/angle_harmonic.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/bonded_coulomb.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/bonded_coulomb_sr.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/bonded_interaction_data.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/bonded_tab.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/dihedral.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/fene.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/harmonic.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/harmonic_dumbbell.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/quartic.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/thermalized_bond.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/umbrella.cpp)
