#ifndef MATHEVAL_IMPLEMENTATION
#error "Do not include ast_adapted.hpp directly!"
#endif

#pragma once

#include "ast.hpp"

#include <boost/fusion/include/adapt_struct.hpp>

BOOST_FUSION_ADAPT_STRUCT(matheval::ast::unary_op, op, rhs)

BOOST_FUSION_ADAPT_STRUCT(matheval::ast::binary_op, op, lhs, rhs)

BOOST_FUSION_ADAPT_STRUCT(matheval::ast::operation, op, rhs)

BOOST_FUSION_ADAPT_STRUCT(matheval::ast::expression, lhs, rhs)
