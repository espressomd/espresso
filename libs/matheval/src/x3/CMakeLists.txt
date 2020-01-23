add_library(matheval.x3 SHARED
  ast_adapted.hpp
  ast.hpp
  evaluator.cpp
  evaluator.hpp
  matheval.cpp
  math.hpp
  parser.cpp
  parser_def.hpp
  parser.hpp
  )
target_include_directories(matheval.x3
  PUBLIC ../../include/matheval/x3
  )
target_compile_features(matheval.x3
  PUBLIC cxx_std_14
  )
find_package(Boost 1.65.1 REQUIRED)
target_link_libraries(matheval.x3
  PRIVATE Boost::boost
  )
add_library(matheval::x3 ALIAS matheval.x3)
