#include <iostream>

#define BOOST_TEST_MODULE AutoParameters test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "script_interface/auto_parameters/AutoParameters.hpp"

using ScriptInterface::AutoParameters;

struct A : AutoParameters {
  A() : AutoParameters({{"i", i}, {"j", j}}), j(42) {}

  const std::string name() const override { return "A"; }

  int i;
  const int j;
  A *p;
};

BOOST_AUTO_TEST_CASE(bla) {
  A a;
  a.i = 0;

  a.set_parameter("i", 5);

  BOOST_CHECK(a.i == boost::get<int>(a.get_parameter("i")));
  BOOST_CHECK(a.j == boost::get<int>(a.get_parameter("j")));
}
