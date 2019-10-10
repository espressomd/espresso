#ifndef MATHEVAL_IMPLEMENTATION
#error "Do not include parser_def.hpp directly!"
#endif

#pragma once

#include "ast.hpp"
#include "ast_adapted.hpp"
#include "math.hpp"
#include "parser.hpp"

#include <boost/math/constants/constants.hpp>
#include <boost/spirit/home/x3.hpp>

#include <cmath>
#include <iostream>
#include <limits>
#include <string>

namespace matheval {

namespace x3 = boost::spirit::x3;

namespace parser {

// LOOKUP

struct constant_ : x3::symbols<double> {
    constant_() {
        add("e", boost::math::constants::e<double>())(
            "epsilon", std::numeric_limits<double>::epsilon())(
            "phi", boost::math::constants::phi<double>())(
            "pi", boost::math::constants::pi<double>());
    }
} constant;

struct ufunc_ : x3::symbols<double (*)(double)> {
    ufunc_() {
        add("abs", static_cast<double (*)(double)>(&std::abs))(
            "acos", static_cast<double (*)(double)>(&std::acos))(
            "acosh", static_cast<double (*)(double)>(&std::acosh))(
            "asin", static_cast<double (*)(double)>(&std::asin))(
            "asinh", static_cast<double (*)(double)>(&std::asinh))(
            "atan", static_cast<double (*)(double)>(&std::atan))(
            "atanh", static_cast<double (*)(double)>(&std::atanh))(
            "cbrt", static_cast<double (*)(double)>(&std::cbrt))(
            "ceil", static_cast<double (*)(double)>(&std::ceil))(
            "cos", static_cast<double (*)(double)>(&std::cos))(
            "cosh", static_cast<double (*)(double)>(&std::cosh))(
            "deg", static_cast<double (*)(double)>(&math::deg))(
            "erf", static_cast<double (*)(double)>(&std::erf))(
            "erfc", static_cast<double (*)(double)>(&std::erfc))(
            "exp", static_cast<double (*)(double)>(&std::exp))(
            "exp2", static_cast<double (*)(double)>(&std::exp2))(
            "floor", static_cast<double (*)(double)>(&std::floor))(
            "isinf", static_cast<double (*)(double)>(&math::isinf))(
            "isnan", static_cast<double (*)(double)>(&math::isnan))(
            "log", static_cast<double (*)(double)>(&std::log))(
            "log2", static_cast<double (*)(double)>(&std::log2))(
            "log10", static_cast<double (*)(double)>(&std::log10))(
            "rad", static_cast<double (*)(double)>(&math::rad))(
            "round", static_cast<double (*)(double)>(&std::round))(
            "sgn", static_cast<double (*)(double)>(&math::sgn))(
            "sin", static_cast<double (*)(double)>(&std::sin))(
            "sinh", static_cast<double (*)(double)>(&std::sinh))(
            "sqrt", static_cast<double (*)(double)>(&std::sqrt))(
            "tan", static_cast<double (*)(double)>(&std::tan))(
            "tanh", static_cast<double (*)(double)>(&std::tanh))(
            "tgamma", static_cast<double (*)(double)>(&std::tgamma));
    }
} ufunc;

struct bfunc_ : x3::symbols<double (*)(double, double)> {
    bfunc_() {
        add("atan2", static_cast<double (*)(double, double)>(&std::atan2))(
            "max", static_cast<double (*)(double, double)>(&std::fmax))(
            "min", static_cast<double (*)(double, double)>(&std::fmin))(
            "pow", static_cast<double (*)(double, double)>(&std::pow));
    }
} bfunc;

struct plusminus_ : x3::symbols<double (*)(double)> {
    plusminus_() {
        add("+", static_cast<double (*)(double)>(&math::plus))(
            "-", static_cast<double (*)(double)>(&math::minus));
    }
} plusminus;

struct addsub_ : x3::symbols<double (*)(double, double)> {
    addsub_() {
        add("+", static_cast<double (*)(double, double)>(&math::plus))(
            "-", static_cast<double (*)(double, double)>(&math::minus));
    }
} addsub;

struct muldiv_ : x3::symbols<double (*)(double, double)> {
    muldiv_() {
        add("*", static_cast<double (*)(double, double)>(&math::multiplies))(
            "/", static_cast<double (*)(double, double)>(&math::divides))(
            "%", static_cast<double (*)(double, double)>(&std::fmod));
    }
} muldiv;

struct power_ : x3::symbols<double (*)(double, double)> {
    power_() { add("**", static_cast<double (*)(double, double)>(&std::pow)); }
} power;

// ADL markers

struct expression_class;
struct term_class;
struct factor_class;
struct primary_class;
struct unary_class;
struct binary_class;
struct variable_class;

// Rule declarations

auto const expression =
    x3::rule<expression_class, ast::expression>{"expression"};
auto const term = x3::rule<term_class, ast::expression>{"term"};
auto const factor = x3::rule<factor_class, ast::expression>{"factor"};
auto const primary = x3::rule<primary_class, ast::operand>{"primary"};
auto const unary = x3::rule<unary_class, ast::unary_op>{"unary"};
auto const binary = x3::rule<binary_class, ast::binary_op>{"binary"};
auto const variable = x3::rule<variable_class, std::string>{"variable"};

// Rule defintions

// clang-format off

auto const expression_def =
    term > *(addsub > term)
    ;

auto const term_def =
    factor > *(muldiv > factor)
    ;

auto const factor_def =
    primary > *( power > factor )
    ;

auto const unary_def =
    ufunc > '(' > expression > ')'
    ;

auto const binary_def =
    bfunc > '(' > expression > ',' > expression > ')'
    ;

auto const variable_def =
    x3::raw[x3::lexeme[x3::alpha >> *(x3::alnum | '_')]]
    ;

auto const primary_def =
      x3::double_
    | ('(' > expression > ')')
    | (plusminus > primary)
    | binary
    | unary
    | constant
    | variable
    ;

// clang-format on

BOOST_SPIRIT_DEFINE(expression, term, factor, primary, unary, binary, variable)

struct expression_class {
    template <typename Iterator, typename Exception, typename Context>
    x3::error_handler_result on_error(Iterator &, Iterator const &last,
                                      Exception const &x, Context const &) {
        std::cout << "Expected " << x.which() << " at \""
                  << std::string{x.where(), last} << "\"" << std::endl;
        return x3::error_handler_result::fail;
    }
};

} // namespace parser

parser::expression_type grammar() { return parser::expression; }

} // namespace matheval
