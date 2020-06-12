# TODO: License header
project(StokesianDynamics LANGUAGES CXX)

add_library(stokesian_dynamics INTERFACE)
target_include_directories(stokesian_dynamics INTERFACE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
add_subdirectory(src)

export(EXPORT stokesiandynamics-targets
       FILE ${CMAKE_CURRENT_BINARY_DIR}/StokesianDynamicsTargets.cmake
       NAMESPACE StokesianDynamics::)

export(PACKAGE StokesianDynamics)
