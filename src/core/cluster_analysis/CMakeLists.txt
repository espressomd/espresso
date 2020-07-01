set(cluster_analysis_SRC Cluster.cpp ClusterStructure.cpp)
add_library(core_cluster_analysis SHARED ${cluster_analysis_SRC})
install(TARGETS core_cluster_analysis
        LIBRARY DESTINATION ${PYTHON_INSTDIR}/espressomd)
set_target_properties(core_cluster_analysis PROPERTIES MACOSX_RPATH TRUE)
target_link_libraries(core_cluster_analysis PUBLIC EspressoCore
                      PRIVATE EspressoConfig cxx_interface)

if(GSL)
  target_link_libraries(core_cluster_analysis PRIVATE GSL::gsl GSL::gslcblas)
endif()
