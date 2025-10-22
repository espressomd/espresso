/*
 * Copyright (C) 2022 The ESPResSo project
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

#define BOOST_TEST_MODULE "Bond breakage"
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "bond_breakage/actions.hpp"

BOOST_AUTO_TEST_CASE(test_actions_equality) {
  {
    using Action = BondBreakage::DeleteBond;
    BOOST_CHECK((Action{1, 2, 3} == Action{1, 2, 3}));
    BOOST_CHECK((Action{1, 2, 3} != Action{0, 2, 3}));
    BOOST_CHECK((Action{1, 2, 3} != Action{1, 0, 3}));
    BOOST_CHECK((Action{1, 2, 3} != Action{1, 2, 0}));
  }

  {
    using Action = BondBreakage::DeleteAngleBond;
    BOOST_CHECK((Action{1, {2, 4}, 3} == Action{1, {2, 4}, 3}));
    BOOST_CHECK((Action{1, {2, 4}, 3} != Action{0, {2, 4}, 3}));
    BOOST_CHECK((Action{1, {2, 4}, 3} != Action{1, {0, 4}, 3}));
    BOOST_CHECK((Action{1, {2, 4}, 3} != Action{1, {2, 0}, 3}));
    BOOST_CHECK((Action{1, {2, 4}, 3} != Action{1, {2, 4}, 0}));
  }

  {
    using Action = BondBreakage::DeleteAllBonds;
    BOOST_CHECK((Action{1, 2} == Action{1, 2}));
    BOOST_CHECK((Action{1, 2} != Action{0, 2}));
    BOOST_CHECK((Action{1, 2} != Action{1, 0}));
  }
}

BOOST_AUTO_TEST_CASE(test_actions_hash_value) {
  {
    using Action = BondBreakage::DeleteBond;
    BOOST_CHECK((Action{1, 2, 3}.hash_value() == Action{1, 2, 3}.hash_value()));
    BOOST_CHECK((Action{1, 2, 3}.hash_value() != Action{0, 2, 3}.hash_value()));
    BOOST_CHECK((Action{1, 2, 3}.hash_value() != Action{1, 0, 3}.hash_value()));
    BOOST_CHECK((Action{1, 2, 3}.hash_value() != Action{1, 2, 0}.hash_value()));
  }

  {
    // clang-format off
    using Action = BondBreakage::DeleteAngleBond;
    BOOST_CHECK((Action{1, {2, 4}, 3}.hash_value() == Action{1, {2, 4}, 3}.hash_value()));
    BOOST_CHECK((Action{1, {2, 4}, 3}.hash_value() != Action{0, {2, 4}, 3}.hash_value()));
    BOOST_CHECK((Action{1, {2, 4}, 3}.hash_value() != Action{1, {0, 4}, 3}.hash_value()));
    BOOST_CHECK((Action{1, {2, 4}, 3}.hash_value() != Action{1, {2, 0}, 3}.hash_value()));
    BOOST_CHECK((Action{1, {2, 4}, 3}.hash_value() != Action{1, {2, 4}, 0}.hash_value()));
    // clang-format on
  }

  {
    using Action = BondBreakage::DeleteAllBonds;
    BOOST_CHECK((Action{1, 2}.hash_value() == Action{1, 2}.hash_value()));
    BOOST_CHECK((Action{1, 2}.hash_value() != Action{0, 2}.hash_value()));
    BOOST_CHECK((Action{1, 2}.hash_value() != Action{1, 0}.hash_value()));
  }
}
