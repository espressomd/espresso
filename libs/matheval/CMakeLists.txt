cmake_minimum_required(VERSION 3.1)
project(matheval LANGUAGES CXX)
find_package(Boost REQUIRED)

add_library(matheval evaluator.cpp matheval.cpp parser.cpp)
set_target_properties(matheval PROPERTIES
  POSITION_INDEPENDENT_CODE ON
  CXX_CLANG_TIDY "")
target_include_directories(matheval PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
target_link_libraries(matheval PUBLIC Boost::boost)
