#define BOOST_TEST_MODULE observables_test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <particle/particle.hpp>
#include <observables/center_of_mass.hpp>

BOOST_AUTO_TEST_CASE(observables) {
  Particle::Particle p1{1.5 /* charge */, 2.0 /* position */, 4.2 /* velocity */, 5.1 /* mass */};
  Particle::Particle p2{1.5, 3.0, 4.2, 3.1};
  std::vector<Particle::Particle> particles{p1, p2};
  {
    auto const res = Observables::CenterOfMass()(particles);;
    BOOST_CHECK(res == (p1.x * p1.m + p2.x * p2.m) / (p1.m + p2.m));
  }
}
