#define BOOST_TEST_MODULE interactions lennard_jones
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "interactions/lennard_jones.hpp"
#include "utils/Vector.hpp"

constexpr auto tol = std::numeric_limits<double>::epsilon();

BOOST_AUTO_TEST_CASE(lennard_jones_test) {
  Interactions::LennardJones<double> lj{1.5, 2.0};
  Utils::Vector3d r12{1, 0, 0};
  BOOST_CHECK_SMALL(lj.energy(2.0), tol);
  BOOST_CHECK_SMALL(lj.energy(std::pow(2.0, 1. / 6.) * 2.0) + 1.5, tol);
  BOOST_CHECK(lj.force(std::pow(2.0, 1. / 6.) * 2.0, r12) ==
              Utils::Vector3d::broadcast(0.0));
}