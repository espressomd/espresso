#define BOOST_TEST_MODULE interactions central_potential
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "interactions/central_potential.hpp"
#include "utils/Vector.hpp"

constexpr auto tol = std::numeric_limits<double>::epsilon();

namespace Testing {

template <class T> struct Linear {
  using value_type = T;

  value_type m_c1;
  value_type m_c2;
  value_type energy(value_type r) const { return -m_c1 * r + m_c2; }
  Utils::Vector<value_type, 3>
  force(value_type r, Utils::Vector<value_type, 3> const &r12) const {
    return r12 * m_c1;
  }
};

} // namespace Testing

BOOST_AUTO_TEST_CASE(central_potential_test) {
  Testing::Linear<double> lin{3.5, 1.0};
  BOOST_CHECK_SMALL(lin.energy(0.0) - 1, tol);
  BOOST_CHECK_SMALL(lin.energy(1.0) - (-3.5 + 1.0), tol);

  constexpr double cutoff = 3.0;
  constexpr double shift = 1.0;
  constexpr double offset = 1.0;

  Interactions::CentralPotential<Testing::Linear<double>> cp_lin(cutoff, offset,
                                                                 shift, lin);
  BOOST_CHECK_SMALL(cp_lin.energy(1.0) - shift - 1.0, tol);
  Utils::Vector3d r12{1, 0, 0};
  BOOST_CHECK(cp_lin.force(0.0, r12) == Utils::Vector3d{});
  BOOST_CHECK(cp_lin.force(1.0, r12) == r12 * 3.5);
}