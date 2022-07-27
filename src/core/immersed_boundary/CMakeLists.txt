target_sources(
  Espresso_core
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/ibm_tribend.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/ibm_triel.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/ibm_volcons.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/ibm_common.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/ImmersedBoundaries.cpp)
