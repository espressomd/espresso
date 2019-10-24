#ifndef MATHEVAL_IMPLEMENTATION
#error "Do not include parser.hpp directly!"
#endif

#pragma once

#include "ast.hpp"

#include <boost/spirit/include/qi.hpp>

#include <iostream>

namespace matheval {

namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;

namespace parser {

struct expectation_handler {
    template <typename>
    struct result {
        typedef void type;
    };

    template <typename Iterator>
    void operator()(Iterator first, Iterator last,
                    boost::spirit::info const &info) const {
        std::stringstream msg;
        msg << "Expected " << info << " at \"" << std::string(first, last)
            << "\"";

        throw std::runtime_error(msg.str()); // NOLINT
    }
};

template <typename Iterator>
struct grammar : qi::grammar<Iterator, ast::expression(), ascii::space_type> {
    expectation_handler err_handler;
    qi::rule<Iterator, ast::expression(), ascii::space_type> expression;
    qi::rule<Iterator, ast::expression(), ascii::space_type> term;
    qi::rule<Iterator, ast::expression(), ascii::space_type> factor;
    qi::rule<Iterator, ast::operand(), ascii::space_type> primary;
    qi::rule<Iterator, ast::unary_op(), ascii::space_type> unary;
    qi::rule<Iterator, ast::binary_op(), ascii::space_type> binary;
    qi::rule<Iterator, std::string()> variable;

    qi::symbols<typename std::iterator_traits<Iterator>::value_type, double>
        constant;
    qi::symbols<typename std::iterator_traits<Iterator>::value_type,
                double (*)(double)>
        ufunc;
    qi::symbols<typename std::iterator_traits<Iterator>::value_type,
                double (*)(double, double)>
        bfunc;
    qi::symbols<typename std::iterator_traits<Iterator>::value_type,
                double (*)(double)>
        plusminus;
    qi::symbols<typename std::iterator_traits<Iterator>::value_type,
                double (*)(double, double)>
        addsub;
    qi::symbols<typename std::iterator_traits<Iterator>::value_type,
                double (*)(double, double)>
        muldiv;
    qi::symbols<typename std::iterator_traits<Iterator>::value_type,
                double (*)(double, double)>
        power;

    grammar();
};

} // namespace parser

typedef parser::grammar<std::string::const_iterator> grammar;

} // namespace matheval
