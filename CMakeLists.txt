cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

project(boost_matheval CXX)
set(CMAKE_CXX_STANDARD 11)

if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the build type" FORCE)
endif(NOT CMAKE_BUILD_TYPE)

enable_testing()
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -C ${CMAKE_BUILD_TYPE})

find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
include_directories(include)

add_subdirectory(doc)
add_subdirectory(examples)
add_subdirectory(tests)

include(FeatureSummary)
feature_summary(WHAT ALL)
