cmake_minimum_required(VERSION 3.10)
project(observables LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 14)
enable_testing()

option(BUILD_TESTING "Build tests" ON)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

add_library(observables INTERFACE)
target_include_directories(observables INTERFACE ${CMAKE_SOURCE_DIR}/include/)

if(BUILD_TESTING)
    add_subdirectory(tests)
endif()
