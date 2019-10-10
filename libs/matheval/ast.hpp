#ifndef MATHEVAL_IMPLEMENTATION
#error "Do not include ast.hpp directly!"
#endif

#pragma once

#include <boost/spirit/home/x3.hpp>
#include <boost/spirit/home/x3/support/ast/variant.hpp>

#include <list>
#include <string>

namespace matheval {

namespace x3 = boost::spirit::x3;

namespace ast {

struct nil {};
struct unary_op;
struct binary_op;
struct expression;

// clang-format off
struct operand : x3::variant<
                 nil
                 , double
                 , std::string
                 , x3::forward_ast<unary_op>
                 , x3::forward_ast<binary_op>
                 , x3::forward_ast<expression>
                 > {
    using base_type::base_type;
    using base_type::operator=;
};
// clang-format on

struct unary_op {
    double (*op)(double);
    operand rhs;
};

struct binary_op {
    double (*op)(double, double);
    operand lhs;
    operand rhs;
};

struct operation {
    double (*op)(double, double);
    operand rhs;
};

struct expression {
    operand lhs;
    std::list<operation> rhs;
};

} // namespace ast

} // namespace matheval
