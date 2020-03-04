#define BOOST_TEST_MODULE traits_test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <fluid/fluid.hpp>
#include <particle/particle.hpp>
#include <traits/fluid.hpp>
#include <traits/particle.hpp>

BOOST_AUTO_TEST_CASE(particle_traits) {
  Particle::Particle p{1.5, 2.0, 4.2, 5.1};
  {
    auto const q = Traits::Particle::Charge()(p);
    BOOST_CHECK(q == 1.5);
  }
  {
    auto const x = Traits::Particle::Position()(p);
    BOOST_CHECK(x == 2.0);
  }
  {
    auto const v = Traits::Particle::Velocity()(p);
    BOOST_CHECK(v == 4.2);
  }
  {
    auto const mass = Traits::Particle::Mass()(p);
    BOOST_CHECK(mass == 5.1);
  }
  {
    auto p2 = p;
    p2.is_virtual = true;
    auto const mass = Traits::Particle::Mass()(p2);
    BOOST_CHECK(mass == 0.0);
  }
}

BOOST_AUTO_TEST_CASE(fluid_traits) {
  Fluid::FluidNode f{2.4};
  auto const res = Traits::Fluid::Velocity()(f);
  BOOST_CHECK(res == 2.4);
}