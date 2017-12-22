#define BOOST_TEST_MODULE AutoParameters test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "script_interface/auto_parameters/AutoParameters.hpp"

using ScriptInterface::AutoParameters;

struct A : AutoParameters {
  A(int i_, int j_) : AutoParameters({{"i", i}, {"j", j}}), i(i_), j(j_) {}

  int i;
  const int j;
};

BOOST_AUTO_TEST_CASE(basic) {
  A a{0, 42};

  auto valid_parameters = a.valid_parameters();
  BOOST_CHECK(valid_parameters.size() == 2);
  BOOST_CHECK(valid_parameters.find("i") != valid_parameters.end());
  BOOST_CHECK(valid_parameters.find("j") != valid_parameters.end());

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
                        AutoParameters::UnknownParameter,
                        [](std::runtime_error const &) { return true; });
  BOOST_CHECK_EXCEPTION(a.set_parameter("unknown", 12),
                        AutoParameters::UnknownParameter,
                        [](std::runtime_error const &) { return true; });
  BOOST_CHECK_EXCEPTION(a.set_parameter("j", 12), AutoParameters::WriteError,
                        [](std::runtime_error const &) { return true; });
}
