set(cluster_analysis_SRC
    Cluster.cpp
    ClusterStructure.cpp
)
add_library(cluster_analysis SHARED ${cluster_analysis_SRC})
install(TARGETS cluster_analysis LIBRARY DESTINATION ${PYTHON_INSTDIR}/espressomd)
set_target_properties(cluster_analysis PROPERTIES SOVERSION ${SOVERSION})
set_target_properties(cluster_analysis PROPERTIES MACOSX_RPATH TRUE)
target_link_libraries(cluster_analysis PUBLIC EspressoCore PRIVATE EspressoConfig)

if(GSL)
target_link_libraries(cluster_analysis PRIVATE GSL::gsl GSL::gslcblas)
endif()
