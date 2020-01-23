#ifndef MATHEVAL_IMPLEMENTATION
#error "Do not include parser_def.hpp directly!"
#endif

#pragma once

#include "ast.hpp"
#include "ast_adapted.hpp"
#include "math.hpp"
#include "parser.hpp"

#include <boost/spirit/include/phoenix.hpp>
#define BOOST_SPIRIT_NO_PREDEFINED_TERMINALS
#include <boost/spirit/include/qi.hpp>

#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

namespace matheval {

namespace qi = boost::spirit::qi;

namespace parser {

template <typename Iterator>
grammar<Iterator>::grammar() : grammar::base_type(expression) {
    qi::_2_type _2;
    qi::_3_type _3;
    qi::_4_type _4;

    qi::alnum_type alnum;
    qi::alpha_type alpha;
    qi::double_type double_;
    qi::lexeme_type lexeme;
    qi::raw_type raw;

    // clang-format off

    constant.add
        ("e"      , boost::math::constants::e<double>())
        ("epsilon", std::numeric_limits<double>::epsilon())
        ("phi"    , boost::math::constants::phi<double>())
        ("pi"     , boost::math::constants::pi<double>())
        ;

    ufunc.add
        ("abs"  , static_cast<double (*)(double)>(&std::abs))
        ("acos" , static_cast<double (*)(double)>(&std::acos))
        ("asin" , static_cast<double (*)(double)>(&std::asin))
        ("atan" , static_cast<double (*)(double)>(&std::atan))
        ("ceil" , static_cast<double (*)(double)>(&std::ceil))
        ("cos"  , static_cast<double (*)(double)>(&std::cos))
        ("cosh" , static_cast<double (*)(double)>(&std::cosh))
        ("deg"  , static_cast<double (*)(double)>(&math::deg))
        ("exp"  , static_cast<double (*)(double)>(&std::exp))
        ("floor", static_cast<double (*)(double)>(&std::floor))
        ("isinf", static_cast<double (*)(double)>(&math::isinf))
        ("isnan", static_cast<double (*)(double)>(&math::isnan))
        ("log"  , static_cast<double (*)(double)>(&std::log))
        ("log10", static_cast<double (*)(double)>(&std::log10))
        ("rad"  , static_cast<double (*)(double)>(&math::rad))
        ("sgn"  , static_cast<double (*)(double)>(&math::sgn))
        ("sin"  , static_cast<double (*)(double)>(&std::sin))
        ("sinh" , static_cast<double (*)(double)>(&std::sinh))
        ("sqrt" , static_cast<double (*)(double)>(&std::sqrt))
        ("tan"  , static_cast<double (*)(double)>(&std::tan))
        ("tanh" , static_cast<double (*)(double)>(&std::tanh))
        ;

    bfunc.add
        ("atan2", static_cast<double (*)(double, double)>(&std::atan2))
        ("pow"  , static_cast<double (*)(double, double)>(&std::pow))
        ;

    unary_op.add
        ("+", static_cast<double (*)(double)>(&math::plus))
        ("-", static_cast<double (*)(double)>(&math::minus))
        ("!", static_cast<double (*)(double)>(&math::unary_not))
        ;

    additive_op.add
        ("+", static_cast<double (*)(double, double)>(&math::plus))
        ("-", static_cast<double (*)(double, double)>(&math::minus))
        ;

    multiplicative_op.add
        ("*", static_cast<double (*)(double, double)>(&math::multiplies))
        ("/", static_cast<double (*)(double, double)>(&math::divides))
        ("%", static_cast<double (*)(double, double)>(&std::fmod))
        ;

    logical_op.add
        ("&&", static_cast<double (*)(double, double)>(&math::logical_and))
        ("||", static_cast<double (*)(double, double)>(&math::logical_or))
        ;

    relational_op.add
        ("<" , static_cast<double (*)(double, double)>(&math::less))
        ("<=", static_cast<double (*)(double, double)>(&math::less_equals))
        (">" , static_cast<double (*)(double, double)>(&math::greater))
        (">=", static_cast<double (*)(double, double)>(&math::greater_equals))
        ;

    equality_op.add
        ("==", static_cast<double (*)(double, double)>(&math::equals))
        ("!=", static_cast<double (*)(double, double)>(&math::not_equals))
        ;

    power.add
        ("**", static_cast<double (*)(double, double)>(&std::pow))
        ;

    expression =
        logical.alias()
        ;

    logical =
        equality >> *(logical_op > equality)
        ;

    equality =
        relational >> *(equality_op > relational)
        ;

    relational =
        additive >> *(relational_op > additive)
        ;

    additive =
        multiplicative >> *(additive_op > multiplicative)
        ;

    multiplicative =
        factor >> *(multiplicative_op > factor)
        ;

    factor =
        primary >> *( power > factor )
        ;

    unary =
        ufunc > '(' > expression > ')'
        ;

    binary =
        bfunc > '(' > expression > ',' > expression > ')'
        ;

    variable =
        raw[lexeme[alpha >> *(alnum | '_')]]
        ;

    primary =
          double_
        | ('(' > expression > ')')
        | (unary_op > primary)
        | binary
        | unary
        | constant
        | variable
        ;

    // clang-format on

    expression.name("expression");
    logical.name("logical");
    equality.name("equality");
    relational.name("relational");
    additive.name("additive");
    multiplicative.name("multiplicative");
    factor.name("factor");
    variable.name("variable");
    primary.name("primary");
    unary.name("unary");
    binary.name("binary");

    // typedef boost::phoenix::function<error_handler<Iterator> >
    // error_handler_function; qi::on_error<qi::fail>(expression,
    //        error_handler_function(error_handler<Iterator>())(
    //            "Error! Expecting ", qi::_4, qi::_3));
    qi::on_error<qi::fail>(
        expression,
        boost::phoenix::bind(boost::phoenix::ref(err_handler), _3, _2, _4));
}

} // namespace parser

} // namespace matheval
