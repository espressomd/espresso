/*
 * Copyright (C) 2026 The ESPResSo project
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#define BOOST_TEST_MODULE TabulatedPotential
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "TabulatedPotential.hpp"

#include <cmath>
#include <vector>

/* A single-point table with minval == maxval is an intentionally supported
 * configuration: the constructor throws only if the force table does not hold
 * exactly one element when minval == maxval (mirrored by the Python test
 * testsuite/python/tabulated.py). With minval == maxval the step size is
 * degenerate (0/0), which previously stored invstepsize = NaN and poisoned
 * every force/energy evaluation while reading out of bounds on the size-1
 * table. The interpolation must instead collapse to the constant table value.
 */
BOOST_AUTO_TEST_CASE(degenerate_single_point_table) {
  // ctor argument order is (min, max, force, energy)
  auto const pot = TabulatedPotential(1., 1., {2.}, {3.});

  // the inverse step size must be a finite number (was 0/0 == NaN). This
  // check holds regardless of whether assertions are enabled.
  BOOST_CHECK(std::isfinite(pot.invstepsize));
  BOOST_CHECK_EQUAL(pot.invstepsize, 0.);

  // evaluating at the single sampling point must return the tabulated value
  // (a one-point table is a constant) and never NaN / an out-of-bounds read.
  BOOST_CHECK(std::isfinite(pot.force(1.)));
  BOOST_CHECK(std::isfinite(pot.energy(1.)));
  BOOST_CHECK_EQUAL(pot.force(1.), 2.);
  BOOST_CHECK_EQUAL(pot.energy(1.), 3.);

  // values outside the range are clamped onto the single point and must also
  // stay finite and constant.
  BOOST_CHECK(std::isfinite(pot.force(0.5)));
  BOOST_CHECK(std::isfinite(pot.energy(2.0)));
  BOOST_CHECK_EQUAL(pot.force(0.5), 2.);
  BOOST_CHECK_EQUAL(pot.energy(2.0), 3.);
}

/* A regular multi-point table must keep behaving correctly. */
BOOST_AUTO_TEST_CASE(regular_table) {
  auto const pot = TabulatedPotential(1., 3., {1., 2., 3.}, {4., 5., 6.});
  BOOST_CHECK(std::isfinite(pot.invstepsize));
  BOOST_CHECK_EQUAL(pot.invstepsize, 1.); // (3-1)/(3-1)
  BOOST_CHECK_EQUAL(pot.force(1.), 1.);
  BOOST_CHECK_EQUAL(pot.force(2.), 2.);
  BOOST_CHECK_CLOSE(pot.force(2.5), 2.5, 1e-12);
  BOOST_CHECK_EQUAL(pot.energy(3.), 6.);
  BOOST_CHECK_EQUAL(pot.cutoff(), 3.);
}
