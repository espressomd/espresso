target_sources(
  EspressoCore PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/Cluster.cpp
                       ${CMAKE_CURRENT_SOURCE_DIR}/ClusterStructure.cpp)

if(GSL)
  target_link_libraries(EspressoCore PRIVATE GSL::gsl GSL::gslcblas)
endif()
