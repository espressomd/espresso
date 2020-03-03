cmake_minimum_required(VERSION 3.10)
project(observables LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 14)
enable_testing()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")
find_package(Boost 1.58 COMPONENTS unit_test_framework REQUIRED)

add_subdirectory(tests)

add_library(traits INTERFACE)
target_include_directories(traits INTERFACE ${CMAKE_SOURCE_DIR}/include/)
add_library(particle INTERFACE)
target_include_directories(particle INTERFACE ${CMAKE_SOURCE_DIR}/include/)
add_library(algorithms INTERFACE)
target_include_directories(algorithms INTERFACE ${CMAKE_SOURCE_DIR}/include/)
add_library(observables INTERFACE)
target_include_directories(observables INTERFACE ${CMAKE_SOURCE_DIR}/include/)
