target_sources(
  Espresso_core PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/Cluster.cpp
                        ${CMAKE_CURRENT_SOURCE_DIR}/ClusterStructure.cpp)

if(GSL)
  target_link_libraries(Espresso_core PRIVATE GSL::gsl GSL::gslcblas)
endif()
