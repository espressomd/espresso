target_sources(
  espresso_core
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/AtomDecomposition.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/CellStructure.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/HybridDecomposition.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/RegularDecomposition.cpp)
