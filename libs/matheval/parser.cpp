#define MATHEVAL_IMPLEMENTATION

#include "parser_def.hpp"

#include <boost/spirit/home/x3.hpp>

#include <string>

namespace matheval {

namespace x3 = boost::spirit::x3;

namespace parser {

using iterator_type = std::string::const_iterator;
using context_type = x3::phrase_parse_context<x3::ascii::space_type>::type;

BOOST_SPIRIT_INSTANTIATE(expression_type, iterator_type, context_type)

} // namespace parser

} // namespace matheval
