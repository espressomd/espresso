/*
 * Copyright (C) 2010-2019 The ESPResSo project
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#define BOOST_TEST_MODULE AutoParameters test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <boost/range/algorithm/find.hpp>

#include "script_interface/auto_parameters/AutoParameters.hpp"

using ScriptInterface::AutoParameters;

struct A : AutoParameters<A> {
  A(int i_, int j_) : AutoParameters({{"i", i}, {"j", j}}), i(i_), j(j_) {}

  int i;
  const int j;
};

BOOST_AUTO_TEST_CASE(basic) {
  A a{0, 42};

  auto const &valid_parameters = a.valid_parameters();

  BOOST_CHECK(valid_parameters.size() == 2);
  BOOST_CHECK(boost::find(valid_parameters, "i") != valid_parameters.end());
  BOOST_CHECK(boost::find(valid_parameters, "j") != valid_parameters.end());

  BOOST_CHECK(0 == boost::get<int>(a.get_parameter("i")));
  BOOST_CHECK(42 == boost::get<int>(a.get_parameter("j")));

  a.set_parameter("i", 12);

  BOOST_CHECK(12 == boost::get<int>(a.get_parameter("i")));
  BOOST_CHECK(42 == boost::get<int>(a.get_parameter("j")));
}

struct B : public A {
  B(int i_, int j_, int k_) : A(i_, j_), k(k_) { add_parameters({{"k", k}}); }
  int k;
};

BOOST_AUTO_TEST_CASE(add_parameters) {
  B b{1, 2, 3};

  BOOST_CHECK(3 == boost::get<int>(b.get_parameter("k")));
  b.set_parameter("k", 12);
  BOOST_CHECK(12 == boost::get<int>(b.get_parameter("k")));
}

BOOST_AUTO_TEST_CASE(exceptions) {
  A a{0, 42};

  BOOST_CHECK_EXCEPTION(a.get_parameter("unknown"),
                        AutoParameters<A>::UnknownParameter,
                        [](std::runtime_error const &) { return true; });
  BOOST_CHECK_EXCEPTION(a.set_parameter("unknown", 12),
                        AutoParameters<A>::UnknownParameter,
                        [](std::runtime_error const &) { return true; });
  BOOST_CHECK_EXCEPTION(a.set_parameter("j", 12), AutoParameters<A>::WriteError,
                        [](std::runtime_error const &) { return true; });
}
