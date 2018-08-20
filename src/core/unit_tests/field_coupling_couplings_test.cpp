#define BOOST_TEST_MODULE AutoParameter test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "core/field_coupling/couplings/Charge.hpp"
#include "core/field_coupling/couplings/Direct.hpp"
#include "core/field_coupling/couplings/Mass.hpp"
#include "core/field_coupling/couplings/Scaled.hpp"
#include "core/field_coupling/couplings/Viscous.hpp"

using namespace FieldCoupling::Coupling;

BOOST_AUTO_TEST_CASE(charge) {
  static_assert(Charge::is_linear, "");

  struct {
    struct {
      const double q = 1.23;
    } p;
  } p;

  BOOST_CHECK((p.p.q * 2.0) == Charge()(p, 2.0));
}

BOOST_AUTO_TEST_CASE(mass) {
  static_assert(Mass::is_linear, "");

  struct {
    struct {
      const double mass = 1.23;
    } p;
  } p;

  BOOST_CHECK((p.p.mass * 2.0) == Mass()(p, 2.0));
}

BOOST_AUTO_TEST_CASE(direct) {
  static_assert(Direct::is_linear, "");

  BOOST_CHECK(5 == Direct()(0, 5));
}

BOOST_AUTO_TEST_CASE(scaled) {
  static_assert(Scaled::is_linear, "");

  auto const scales = std::unordered_map<int, double>{{0, 1.23}, {2, 3.45}};
  auto const default_val = 9.8;

  auto const scaled_coupling = Scaled(scales, default_val);

  BOOST_CHECK(scales == scaled_coupling.particle_scales());
  BOOST_CHECK(default_val == scaled_coupling.default_scale());

  {
    struct Particle {
      Particle(int id) : m_id(id) {}

      int identity() const { return m_id; };

      const int m_id;
    };

    BOOST_CHECK((1.23 * 2.) == scaled_coupling(Particle(0), 2.));
    BOOST_CHECK((default_val * 3.) == scaled_coupling(Particle(1), 3.));
    BOOST_CHECK((3.45 * 4.) == scaled_coupling(Particle(2), 4.));
    BOOST_CHECK((default_val * 5.) == scaled_coupling(Particle(3), 5.));
  }
}

BOOST_AUTO_TEST_CASE(viscous) {
  static_assert(Viscous::is_linear, "");

  auto const gamma = 3.14159;

  const auto viscous_coupling = Viscous(gamma);

  BOOST_CHECK(gamma == viscous_coupling.gamma());

  {
    struct {
      struct {
        const Vector3d v = {1., 2., 3.};
      } m;
    } p;

    auto const u = Vector3d{4., 5., 6.};

    BOOST_CHECK((-gamma * (p.m.v - u)) == viscous_coupling(p, u));
  }
}
