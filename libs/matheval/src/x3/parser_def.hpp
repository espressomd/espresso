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
        // clang-format off
        add
            ("e"      , boost::math::constants::e<double>())
            ("epsilon", std::numeric_limits<double>::epsilon())
            ("phi"    , boost::math::constants::phi<double>())
            ("pi"     , boost::math::constants::pi<double>())
            ;
        // clang-format on
    }
} constant;

struct ufunc_ : x3::symbols<double (*)(double)> {
    ufunc_() {
        // clang-format off
        add
            ("abs"   , static_cast<double (*)(double)>(&std::abs))
            ("acos"  , static_cast<double (*)(double)>(&std::acos))
            ("acosh" , static_cast<double (*)(double)>(&std::acosh))
            ("asin"  , static_cast<double (*)(double)>(&std::asin))
            ("asinh" , static_cast<double (*)(double)>(&std::asinh))
            ("atan"  , static_cast<double (*)(double)>(&std::atan))
            ("atanh" , static_cast<double (*)(double)>(&std::atanh))
            ("cbrt"  , static_cast<double (*)(double)>(&std::cbrt))
            ("ceil"  , static_cast<double (*)(double)>(&std::ceil))
            ("cos"   , static_cast<double (*)(double)>(&std::cos))
            ("cosh"  , static_cast<double (*)(double)>(&std::cosh))
            ("deg"   , static_cast<double (*)(double)>(&math::deg))
            ("erf"   , static_cast<double (*)(double)>(&std::erf))
            ("erfc"  , static_cast<double (*)(double)>(&std::erfc))
            ("exp"   , static_cast<double (*)(double)>(&std::exp))
            ("exp2"  , static_cast<double (*)(double)>(&std::exp2))
            ("floor" , static_cast<double (*)(double)>(&std::floor))
            ("isinf" , static_cast<double (*)(double)>(&math::isinf))
            ("isnan" , static_cast<double (*)(double)>(&math::isnan))
            ("log"   , static_cast<double (*)(double)>(&std::log))
            ("log2"  , static_cast<double (*)(double)>(&std::log2))
            ("log10" , static_cast<double (*)(double)>(&std::log10))
            ("rad"   , static_cast<double (*)(double)>(&math::rad))
            ("round" , static_cast<double (*)(double)>(&std::round))
            ("sgn"   , static_cast<double (*)(double)>(&math::sgn))
            ("sin"   , static_cast<double (*)(double)>(&std::sin))
            ("sinh"  , static_cast<double (*)(double)>(&std::sinh))
            ("sqrt"  , static_cast<double (*)(double)>(&std::sqrt))
            ("tan"   , static_cast<double (*)(double)>(&std::tan))
            ("tanh"  , static_cast<double (*)(double)>(&std::tanh))
            ("tgamma", static_cast<double (*)(double)>(&std::tgamma))
            ;
        // clang-format on
    }
} ufunc;

struct bfunc_ : x3::symbols<double (*)(double, double)> {
    bfunc_() {
        // clang-format off
        add
            ("atan2", static_cast<double (*)(double, double)>(&std::atan2))
            ("max"  , static_cast<double (*)(double, double)>(&std::fmax))
            ("min"  , static_cast<double (*)(double, double)>(&std::fmin))
            ("pow"  , static_cast<double (*)(double, double)>(&std::pow))
            ;
        // clang-format on
    }
} bfunc;

struct unary_op_ : x3::symbols<double (*)(double)> {
    unary_op_() {
        // clang-format off
        add
            ("+", static_cast<double (*)(double)>(&math::plus))
            ("-", static_cast<double (*)(double)>(&math::minus))
            ("!", static_cast<double (*)(double)>(&math::unary_not))
            ;
        // clang-format on
    }
} unary_op;

struct additive_op_ : x3::symbols<double (*)(double, double)> {
    additive_op_() {
        // clang-format off
        add
            ("+", static_cast<double (*)(double, double)>(&math::plus))
            ("-", static_cast<double (*)(double, double)>(&math::minus))
            ;
        // clang-format on
    }
} additive_op;

struct multiplicative_op_ : x3::symbols<double (*)(double, double)> {
    multiplicative_op_() {
        // clang-format off
        add
            ("*", static_cast<double (*)(double, double)>(&math::multiplies))
            ("/", static_cast<double (*)(double, double)>(&math::divides))
            ("%", static_cast<double (*)(double, double)>(&std::fmod))
            ;
        // clang-format on
    }
} multiplicative_op;

struct logical_op_ : x3::symbols<double (*)(double, double)> {
    logical_op_() {
        // clang-format off
        add
            ("&&", static_cast<double (*)(double, double)>(&math::logical_and))
            ("||", static_cast<double (*)(double, double)>(&math::logical_or))
            ;
        // clang-format on
    }
} logical_op;

struct relational_op_ : x3::symbols<double (*)(double, double)> {
    relational_op_() {
        // clang-format off
        add
            ("<" , static_cast<double (*)(double, double)>(&math::less))
            ("<=", static_cast<double (*)(double, double)>(&math::less_equals))
            (">" , static_cast<double (*)(double, double)>(&math::greater))
            (">=", static_cast<double (*)(double, double)>(&math::greater_equals))
            ;
        // clang-format on
    }
} relational_op;

struct equality_op_ : x3::symbols<double (*)(double, double)> {
    equality_op_() {
        // clang-format off
        add
            ("==", static_cast<double (*)(double, double)>(&math::equals))
            ("!=", static_cast<double (*)(double, double)>(&math::not_equals))
            ;
        // clang-format on
    }
} equality_op;

struct power_ : x3::symbols<double (*)(double, double)> {
    power_() {
        // clang-format off
        add
            ("**", static_cast<double (*)(double, double)>(&std::pow))
            ;
        // clang-format on
    }
} power;

// ADL markers

struct expression_class;
struct logical_class;
struct equality_class;
struct relational_class;
struct additive_class;
struct multiplicative_class;
struct factor_class;
struct primary_class;
struct unary_class;
struct binary_class;
struct variable_class;

// clang-format off

// Rule declarations

auto const expression     = x3::rule<expression_class    , ast::expression>{"expression"};
auto const logical        = x3::rule<logical_class       , ast::expression>{"logical"};
auto const equality       = x3::rule<equality_class      , ast::expression>{"equality"};
auto const relational     = x3::rule<relational_class    , ast::expression>{"relational"};
auto const additive       = x3::rule<additive_class      , ast::expression>{"additive"};
auto const multiplicative = x3::rule<multiplicative_class, ast::expression>{"multiplicative"};
auto const factor         = x3::rule<factor_class        , ast::expression>{"factor"};
auto const primary        = x3::rule<primary_class       , ast::operand   >{"primary"};
auto const unary          = x3::rule<unary_class         , ast::unary_op  >{"unary"};
auto const binary         = x3::rule<binary_class        , ast::binary_op >{"binary"};
auto const variable       = x3::rule<variable_class      , std::string    >{"variable"};

// Rule defintions

auto const expression_def =
    logical
    ;

auto const logical_def =
    equality >> *(logical_op > equality)
    ;

auto const equality_def =
    relational >> *(equality_op > relational)
    ;

auto const relational_def =
    additive >> *(relational_op > additive)
    ;

auto const additive_def =
    multiplicative >> *(additive_op > multiplicative)
    ;

auto const multiplicative_def =
    factor >> *(multiplicative_op > factor)
    ;

auto const factor_def =
    primary >> *( power > factor )
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
    | (unary_op > primary)
    | binary
    | unary
    | constant
    | variable
    ;

BOOST_SPIRIT_DEFINE(
    expression,
    logical,
    equality,
    relational,
    additive,
    multiplicative,
    factor,
    primary,
    unary,
    binary,
    variable
)

// clang-format on

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
