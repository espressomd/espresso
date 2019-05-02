/*
Copyright (C) 2010-2018 The ESPResSo project

This file is part of ESPResSo.

ESPResSo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ESPResSo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
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
