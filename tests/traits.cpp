#define BOOST_TEST_MODULE traits_test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <particle/particle.hpp>
#include <traits/particle.hpp>

BOOST_AUTO_TEST_CASE(traits) {
  Particle::Particle p{1.5, 2.0, 4.2, 5.1};
  {
    auto const res = Traits::Charge()(p);
    BOOST_CHECK(res == 1.5);
  }
  {
    auto const res = Traits::Position()(p);
    BOOST_CHECK(res == 2.0);
  }
  {
    auto const res = Traits::Velocity()(p);
    BOOST_CHECK(res == 4.2);
  }
  {
    auto const res = Traits::Mass()(p);
    BOOST_CHECK(res == 5.1);
  }
}
