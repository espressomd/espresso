#define BOOST_TEST_MODULE observables_test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <observables/center_of_mass.hpp>
#include <observables/current.hpp>
#include <observables/dipole_moment.hpp>
#include <observables/particle.hpp>

BOOST_AUTO_TEST_CASE(observables) {
  Particle::Particle p1{1.5 /* charge */, 2.0 /* position */,
                        4.2 /* velocity */, 5.1 /* mass */};
  Particle::Particle p2{1.5, 3.0, 4.2, 3.1};
  std::vector<Particle::Particle> particles{p1, p2};
  {
    auto const com = Observables::CenterOfMass()(particles);
    BOOST_CHECK(com == (p1.x * p1.m + p2.x * p2.m) / (p1.m + p2.m));
  }
  {
    auto const dipole_moment = Observables::DipoleMoment()(particles);
    BOOST_CHECK(dipole_moment == (p1.x * p1.q + p2.x * p2.q));
  }
  {
    auto const current = Observables::Current()(particles);
    BOOST_CHECK(current == (p1.v * p1.q + p2.v * p2.q));
  }
}
