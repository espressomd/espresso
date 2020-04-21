#define BOOST_TEST_MODULE traits_test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <genobs/observable.hpp>
#include <genobs/properties.hpp>

#include "mock.hpp"

/* Check that the property functors correctly map to
 * the default traits */
BOOST_AUTO_TEST_CASE(properties_) {
  using namespace Observables;

  Testing::Particle p;

  BOOST_CHECK_EQUAL(Position{}(p), traits<Testing::Particle>{}.position(p));
  BOOST_CHECK_EQUAL(Velocity{}(p), traits<Testing::Particle>{}.velocity(p));
  BOOST_CHECK_EQUAL(Mass{}(p), traits<Testing::Particle>{}.mass(p));
  BOOST_CHECK_EQUAL(Charge{}(p), traits<Testing::Particle>{}.charge(p));
  BOOST_CHECK_EQUAL(Force{}(p), traits<Testing::Particle>{}.force(p));
}
