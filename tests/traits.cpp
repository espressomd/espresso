#define BOOST_TEST_MODULE traits_test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <particle/particle.hpp>
#include <fluid/fluid.hpp>
#include <traits/particle.hpp>
#include <traits/fluid.hpp>

BOOST_AUTO_TEST_CASE(particle_traits) {
  Particle::Particle p{1.5, 2.0, 4.2, 5.1};
  {
    auto const res = Traits::Particle::Charge()(p);
    BOOST_CHECK(res == 1.5);
  }
  {
    auto const res = Traits::Particle::Position()(p);
    BOOST_CHECK(res == 2.0);
  }
  {
    auto const res = Traits::Particle::Velocity()(p);
    BOOST_CHECK(res == 4.2);
  }
  {
    auto const res = Traits::Particle::Mass()(p);
    BOOST_CHECK(res == 5.1);
  }
}

BOOST_AUTO_TEST_CASE(fluid_traits) {
  Fluid::FluidNode f{2.4};
  auto const res = Traits::Fluid::Velocity()(f);
  BOOST_CHECK(res == 2.4);
}