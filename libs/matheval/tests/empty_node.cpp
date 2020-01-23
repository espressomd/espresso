#define BOOST_TEST_MODULE empty node
#include <boost/test/included/unit_test.hpp>

#define MATHEVAL_IMPLEMENTATION

#include "ast.hpp"
#include "evaluator.hpp"

BOOST_AUTO_TEST_CASE(empty_node)
{
    matheval::ast::expression ast;
    matheval::ast::eval solver({});
    BOOST_CHECK_EQUAL(solver(ast),0);
}
