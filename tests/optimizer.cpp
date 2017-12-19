#define BOOST_TEST_MODULE optimizer
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "matheval.hpp"

using namespace matheval;

BOOST_AUTO_TEST_CASE(ConstantFolder) {
  /* Binary op */
  {
    auto expr = std::string("1 + 1");
    auto ast = detail::parse<double>(expr.begin(), expr.end());
    auto transformed =
        boost::apply_visitor(detail::ConstantFolder<double>{}, ast.tree);
    BOOST_CHECK_EQUAL(boost::get<double>(transformed), 1. + 1.);
  }

  /* Binary op update right child */
  {
    auto expr = std::string("x + 1 * ( 2 + 3 * 4)");
    auto ast = detail::parse<double>(expr.begin(), expr.end());
    auto transformed =
        boost::apply_visitor(detail::ConstantFolder<double>{}, ast.tree);

    /* x + 1 * ( 2 + 3 * 4) -> x + 14 */
    BOOST_REQUIRE(
        detail::holds_alternative<detail::binary_op<double>>(transformed));
    auto root = boost::get<detail::binary_op<double>>(transformed);
    BOOST_CHECK_EQUAL(1. * (2. + 3. * 4.),
                         boost::get<double>(root.rhs.tree));
  }

  /* Binary op update left child */
  {
    auto expr = std::string("1 * ( 2 + 3 * 4) + x");
    auto ast = detail::parse<double>(expr.begin(), expr.end());
    auto transformed =
      boost::apply_visitor(detail::ConstantFolder<double>{}, ast.tree);

    /* 1 * ( 2 + 3 * 4) + x -> 14 + x */
    BOOST_REQUIRE(
        detail::holds_alternative<detail::binary_op<double>>(transformed));
    auto root = boost::get<detail::binary_op<double>>(transformed);
    BOOST_CHECK_EQUAL(1. * (2. + 3. * 4.), boost::get<double>(root.lhs.tree));
  }

  /* Unary op */
  {
    auto expr = std::string("-(1)");
    auto ast = detail::parse<double>(expr.begin(), expr.end());
    auto transformed =
        boost::apply_visitor(detail::ConstantFolder<double>{}, ast.tree);

    /* -literal(1) -> literal(-1) */
    BOOST_CHECK_EQUAL(-1., boost::get<double>(transformed));
  }

  /* Unary update child */
  {
    auto expr = std::string("-(x + 2 * (1 + 1))");
    auto ast = detail::parse<double>(expr.begin(), expr.end());
    auto transformed =
        boost::apply_visitor(detail::ConstantFolder<double>{}, ast.tree);

    /* -(x + 2 * (1 + 1)) -> -(x + 4) */
    BOOST_REQUIRE(
        detail::holds_alternative<detail::unary_op<double>>(transformed));
    auto minus = boost::get<detail::unary_op<double>>(transformed);
    BOOST_REQUIRE(
        detail::holds_alternative<detail::binary_op<double>>(minus.rhs.tree));
    auto plus = boost::get<detail::binary_op<double>>(minus.rhs.tree);
    BOOST_CHECK_EQUAL("x", boost::get<std::string>(plus.lhs.tree));
    BOOST_CHECK_EQUAL(2. * (1. + 1.), boost::get<double>(plus.rhs.tree));
  }

  /* Variable unchanged */
  {
    auto expr = std::string("x");
    auto ast = detail::parse<double>(expr.begin(), expr.end());
    auto transformed =
        boost::apply_visitor(detail::ConstantFolder<double>{}, ast.tree);

    BOOST_CHECK_EQUAL("x", boost::get<std::string>(transformed));
  }

  /* nil -> 0 */
  {
    detail::expr_ast<double>::tree_t ast = detail::nil{};
    auto transformed =
        boost::apply_visitor(detail::ConstantFolder<double>{}, ast);

    BOOST_CHECK_EQUAL(0, boost::get<double>(transformed));
  }

  /* Evaluate AST */
  {
    auto expr = std::string("1 + 1");
    auto ast = detail::parse<double>(expr.begin(), expr.end());
    auto transformed = detail::ConstantFolder<double>{}(ast);
    BOOST_CHECK_EQUAL(boost::get<double>(transformed), 1. + 1.);
  }
}
