add_library(matheval.x3
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
set_target_properties(matheval.x3
    PROPERTIES CXX_STANDARD 14
  )
add_library(matheval::x3 ALIAS matheval.x3)
