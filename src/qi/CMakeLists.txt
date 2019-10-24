add_library(matheval.qi
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
target_include_directories(matheval.qi
  PUBLIC ../../include/matheval/qi
  )
set_target_properties(matheval.qi
    PROPERTIES CXX_STANDARD 98
  )
add_library(matheval::qi ALIAS matheval.qi)
