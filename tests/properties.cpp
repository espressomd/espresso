#define BOOST_TEST_MODULE traits_test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <observables/properties.hpp>

namespace Testing {
struct Particle {
  double position = 1.;
  double m_vel = 2.;

  auto const &velocity() const { return m_vel; }
};
} // namespace Testing

namespace Observables {
namespace Properties {
template <> struct traits<Testing::Particle> {
  using Particle = Testing::Particle;

  double position(Particle const &p) const { return p.position; }
  double velocity(Particle const &p) const { return p.velocity(); }
  double mass(Particle const &) const { return 1.; }
  double charge(Particle const &) const { return 0.; }
};
} // namespace Properties
} // namespace Observables

/* Check that the Properties::* functors correctly map to
 * the default traits */
BOOST_AUTO_TEST_CASE(properties_) {
  using namespace Observables::Properties;

  Testing::Particle p;

  BOOST_CHECK_EQUAL(Position{}(p), traits<Testing::Particle>{}.position(p));
  BOOST_CHECK_EQUAL(Velocity{}(p), traits<Testing::Particle>{}.velocity(p));
  BOOST_CHECK_EQUAL(Mass{}(p), traits<Testing::Particle>{}.mass(p));
  BOOST_CHECK_EQUAL(Charge{}(p), traits<Testing::Particle>{}.charge(p));
}
