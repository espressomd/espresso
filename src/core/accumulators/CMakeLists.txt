target_sources(
  espresso_core
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/Correlator.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/MeanVarianceCalculator.cpp
          ${CMAKE_CURRENT_SOURCE_DIR}/TimeSeries.cpp)
