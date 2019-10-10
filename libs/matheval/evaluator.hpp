#ifndef MATHEVAL_IMPLEMENTATION
#error "Do not include evaluator.hpp directly!"
#endif

#pragma once

#include "ast.hpp"

#include <map>
#include <string>

namespace matheval {

namespace ast {

struct ConstantFolder {
    using result_type = operand;

    result_type operator()(nil) const;

    result_type operator()(double n) const;

    result_type operator()(std::string const &c) const;

    result_type operator()(operation const &x, operand const &lhs) const;

    result_type operator()(unary_op const &x) const;

    result_type operator()(binary_op const &x) const;

    result_type operator()(expression const &x) const;
};

struct eval {
    using result_type = double;

    explicit eval(std::map<std::string, double> sym) : st(std::move(sym)) {}

    double operator()(nil) const;

    double operator()(double n) const;

    double operator()(std::string const &c) const;

    double operator()(operation const &x, double lhs) const;

    double operator()(unary_op const &x) const;

    double operator()(binary_op const &x) const;

    double operator()(expression const &x) const;

private:
    std::map<std::string, double> st;
};

} // namespace ast

} // namespace matheval
