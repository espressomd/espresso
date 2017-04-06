#include <functional>

#define BOOST_TEST_MODULE make_function test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "utils/make_function.hpp"
using Utils::make_function;

BOOST_AUTO_TEST_CASE(lambda) {
  std::function<int(int)> f = make_function([](int i) { return i + 42; });

  BOOST_CHECK(f(4) == 46);
}

BOOST_AUTO_TEST_CASE(lambda_mutable) {
  std::function<int(int)> f =
      make_function([](int i) mutable { return i + 42; });

  BOOST_CHECK(f(4) == 46);
}

BOOST_AUTO_TEST_CASE(functor) {
  struct F {
    F(int j) : j(j){};
    int operator()(int i) const { return i + j; }

    int j;
  };

  std::function<int(int)> f = make_function(F{42});

  BOOST_CHECK(f(4) == 46);
}

BOOST_AUTO_TEST_CASE(functor_mutable) {
  struct F {
    F(int j) : j(j){};
    int operator()(int i) { return i + j; }

    int j;
  };

  std::function<int(int)> f = make_function(F{42});

  BOOST_CHECK(f(4) == 46);
}

BOOST_AUTO_TEST_CASE(std_function) {
  std::function<int(int)> fun = [](int i) { return 42 + i; };

  std::function<int(int)> f = make_function(fun);

  BOOST_CHECK(f(4) == 46);
}
