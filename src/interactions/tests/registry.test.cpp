#define BOOST_TEST_MODULE interactions registry
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/variant.hpp>

#include <unordered_map>
#include <vector>

#include "interactions/central_potential.hpp"
#include "interactions/lennard_jones.hpp"
#include "interactions/linear.hpp"
#include "interactions/registry.hpp"

#include "utils/Vector.hpp"

BOOST_AUTO_TEST_CASE(registry) {
  using namespace Interactions;
  NonBondedRegistry<Storage<double>, double> registry;
  CentralPotential<Linear<double>> linear(1.0, 0.0, 1.0,
                                          Linear<double>{3.5, 1.0});
  CentralPotential<LennardJones<double>> lj(2.0, 0.0, 1.0,
                                            LennardJones<double>{3.5, 1.0});
  registry.add(0, 0, linear);
  registry.add(0, 0, lj);
  BOOST_CHECK_SMALL(registry.max_cut() - 2.0,
                    std::numeric_limits<double>::epsilon());
}