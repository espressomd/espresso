#define BOOST_TEST_MODULE interactions central_potential
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "interactions/central_potential.hpp"
#include "interactions/linear.hpp"
#include "utils/Vector.hpp"

constexpr auto tol = std::numeric_limits<double>::epsilon();

BOOST_AUTO_TEST_CASE(central_potential_test) {
  Interactions::Linear<double> lin{3.5, 1.0};

  constexpr double cutoff = 3.0;
  constexpr double shift = 1.0;
  constexpr double offset = 1.0;

  Interactions::CentralPotential<Interactions::Linear<double>> cp_lin(
      cutoff, offset, shift, lin);
  BOOST_CHECK_SMALL(cp_lin.energy(1.0) - shift - 1.0, tol);
  Utils::Vector3d r12{1, 0, 0};
  BOOST_CHECK(cp_lin.force(0.0, r12) == Utils::Vector3d{});
  BOOST_CHECK(cp_lin.force(1.0, r12) == r12 * 3.5);
}