#ifndef MATHEVAL_IMPLEMENTATION
#error "Do not include ast_adapted.hpp directly!"
#endif

#pragma once

#include "ast.hpp"

#include <boost/fusion/include/adapt_struct.hpp>

typedef double (*unary_fn)(double);
typedef double (*binary_fn)(double, double);

BOOST_FUSION_ADAPT_STRUCT(matheval::ast::unary_op,
                          (unary_fn, op)(matheval::ast::operand, rhs))

BOOST_FUSION_ADAPT_STRUCT(matheval::ast::binary_op,
                          (binary_fn, op)(matheval::ast::operand,
                                          lhs)(matheval::ast::operand, rhs))

BOOST_FUSION_ADAPT_STRUCT(matheval::ast::operation,
                          (binary_fn, op)(matheval::ast::operand, rhs))

BOOST_FUSION_ADAPT_STRUCT(matheval::ast::expression,
                          (matheval::ast::operand,
                           lhs)(std::list<matheval::ast::operation>, rhs))
