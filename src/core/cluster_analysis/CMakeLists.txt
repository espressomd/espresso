set(cluster_analysis_SRC
    Cluster.cpp
    ClusterStructure.cpp
)
add_library(cluster_analysis SHARED ${cluster_analysis_SRC})
install(TARGETS cluster_analysis LIBRARY DESTINATION ${LIBDIR})
set_target_properties(cluster_analysis PROPERTIES SOVERSION ${SOVERSION})
set_target_properties(cluster_analysis PROPERTIES MACOSX_RPATH TRUE)
add_dependencies(cluster_analysis EspressoConfig)
target_include_directories(cluster_analysis PRIVATE ${CMAKE_SOURCE_DIR}/src/core ${CMAKE_BINARY_DIR}/src/core)
