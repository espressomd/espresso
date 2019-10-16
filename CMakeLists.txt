cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

project(boost_matheval CXX)

# Default build type is Debug
if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the build type" FORCE)
endif()

# Add local package finders
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Build options
option(WITH_COVERAGE "Generate code coverage report" OFF)
if(WITH_COVERAGE)
  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
    link_libraries(gcov)
  else()
    message(WARNING "Coverage information only supported on GCC")
  endif()
endif()

if(CMAKE_VERSION VERSION_GREATER 3.5.2 AND CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  option(WITH_CLANG_TIDY "Run Clang-Tidy during compilation" OFF)
endif()
if(WITH_CLANG_TIDY)
  string(REGEX REPLACE "^([1-9]+\\.[0-9]+).*$" "\\1" CLANG_MINOR_VERSION "${CMAKE_CXX_COMPILER_VERSION}")
  find_program(
    CLANG_TIDY_EXE
    NAMES "clang-tidy-${CLANG_MINOR_VERSION}" "clang-tidy"
    DOC "Path to clang-tidy executable"
    )
  if(NOT CLANG_TIDY_EXE)
    message(STATUS "clang-tidy not found.")
  else()
    message(STATUS "clang-tidy found: ${CLANG_TIDY_EXE}")
    set(CMAKE_CXX_CLANG_TIDY "${CLANG_TIDY_EXE}"
      "-checks=-*,bugprone-*,cert-*,clang-analyzer-*,cppcoreguidelines-*,hicpp-*,misc-*,modernize-*,performance-*,readability-*"
      "-header-filter=.*")
      #"-warnings-as-errors=*")
  endif()
endif()

# Libraries

add_library(matheval_qi
  src/qi/ast_adapted.hpp
  src/qi/ast.hpp
  src/qi/evaluator.cpp
  src/qi/evaluator.hpp
  src/qi/matheval.cpp
  src/qi/math.hpp
  src/qi/parser.cpp
  src/qi/parser_def.hpp
  src/qi/parser.hpp
  )
target_include_directories(matheval_qi
  PUBLIC include/qi
  PRIVATE src/qi
  )
set_target_properties(matheval_qi
  PROPERTIES CXX_STANDARD 98
  )

add_library(matheval_x3
  src/x3/ast_adapted.hpp
  src/x3/ast.hpp
  src/x3/evaluator.cpp
  src/x3/evaluator.hpp
  src/x3/matheval.cpp
  src/x3/math.hpp
  src/x3/parser.cpp
  src/x3/parser_def.hpp
  src/x3/parser.hpp
  )
target_include_directories(matheval_x3
  PUBLIC include/x3
  PRIVATE src/x3
  )
set_target_properties(matheval_x3
  PROPERTIES CXX_STANDARD 14
  )

# Custom test command
enable_testing()
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --output-on-failure -C ${CMAKE_BUILD_TYPE})

include_directories(include/x3)
link_libraries(matheval_x3)

# Find Boost
find_package(Boost REQUIRED COMPONENTS unit_test_framework)
include_directories(${Boost_INCLUDE_DIRS})

add_subdirectory(doc)
add_subdirectory(examples)
add_subdirectory(tests)

include(FeatureSummary)
feature_summary(WHAT ALL)
