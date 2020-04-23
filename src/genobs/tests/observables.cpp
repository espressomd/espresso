#define BOOST_TEST_MODULE observables_test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <vector>

#include <genobs/observable.hpp>

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

  auto prod = GenObs::Product<strip_args<std::integral_constant<int, 2>>,
                              strip_args<std::integral_constant<int, 3>>>{};

  BOOST_CHECK_EQUAL((prod.template operator()<int, int>(0)), 2 * 3);
}

BOOST_AUTO_TEST_CASE(obs) {
  using namespace GenObs;
  Testing::Particle p;
  BOOST_CHECK_EQUAL(Momentum{}(p), Mass{}(p)*Velocity{}(p));
  std::vector<Testing::Particle> parts{p, Testing::Particle{}};
  parts[1].mass = 5.;
  {
    auto const res = AverageMomentum{}(parts);
    BOOST_CHECK(res == 0.5 * (Momentum{}(parts[0]) + Momentum{}(parts[1])));
  }
  {
    auto const res = CenterOfMassPosition{}(parts);
    BOOST_CHECK(res == (Mass{}(parts[0]) * Position{}(parts[0]) +
                        Mass{}(parts[1]) * Position{}(parts[1])) /
                           (Mass{}(parts[0]) + Mass{}(parts[1])));
  }
  {
    auto const res = CenterOfMassVelocity{}(parts);
    BOOST_CHECK(res == (Mass{}(parts[0]) * Velocity{}(parts[0]) +
                        Mass{}(parts[1]) * Velocity{}(parts[1])) /
                           (Mass{}(parts[0]) + Mass{}(parts[1])));
  }
  {
    parts[0].m_force = 1.3;
    parts[1].m_force = 1.8;
    auto const res = TotalForce{}(parts);
    BOOST_CHECK(res == Force{}(parts[0]) + Force{}(parts[1]));
  }
  {
    auto const res = Positions{}(parts);
    BOOST_CHECK(res[0] == Position{}(parts[0]));
    BOOST_CHECK(res[1] == Position{}(parts[1]));
  }
  {
    auto const res = Velocities{}(parts);
    BOOST_CHECK(res[0] == Velocity{}(parts[0]));
    BOOST_CHECK(res[1] == Velocity{}(parts[1]));
  }
  {
    auto const res = Forces{}(parts);
    BOOST_CHECK(res[0] == Force{}(parts[0]));
    BOOST_CHECK(res[1] == Force{}(parts[1]));
  }
}
