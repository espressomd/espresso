target_sources(
  EspressoCore
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/angle_cosine.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/angle_cossquare.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/bonded_interaction_data.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/bonded_tab.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/fene.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/rigid_bond.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/thermalized_bond.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/thermalized_bond_utils.cpp)
