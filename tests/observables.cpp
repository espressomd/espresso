#define BOOST_TEST_MODULE observables_test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <observables/observable.hpp>

#include "mock.hpp"

namespace Testing {
template <class T> struct strip_args {
  template <class... Args> decltype(auto) operator()(Args...) const {
    return T{}();
  }
};
} // namespace Testing

BOOST_AUTO_TEST_CASE(product_) {
  using Testing::strip_args;

  auto prod =
      Observables::Product<strip_args<std::integral_constant<int, 2>>,
                           strip_args<std::integral_constant<int, 3>>>{};

  BOOST_CHECK_EQUAL((prod.template operator()<int, int>(0)), 2 * 3);
}

BOOST_AUTO_TEST_CASE(obs) {
  using namespace Observables;
  Testing::Particle p;
  BOOST_CHECK_EQUAL(Momentum{}(p), Mass{}(p)*Velocity{}(p));
}
