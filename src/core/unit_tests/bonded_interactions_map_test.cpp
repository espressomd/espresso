/*
 * Copyright (C) 2021 The ESPResSo project
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

#define BOOST_TEST_MODULE BondedInteractionsMap test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "bonded_interactions/angle_cosine.hpp"
#include "bonded_interactions/bonded_interaction_data.hpp"
#include "bonded_interactions/fene.hpp"

#include <unordered_map>

BOOST_AUTO_TEST_CASE(insert_bond_types) {
  BondedInteractionsMap bond_map{};
  std::unordered_map<int, std::shared_ptr<Bonded_IA_Parameters>> mock_core{};
  // check defaulted maps are empty
  BOOST_TEST(bond_map.empty());
  BOOST_TEST(mock_core.empty());
  // insert first element
  int first_key = 1;
  auto const fene_bond = FeneBond(1.2, 3.0, 0.8);
  auto const fene_bond_ia = std::make_shared<Bonded_IA_Parameters>(fene_bond);
  bond_map.insert(first_key, fene_bond_ia);
  mock_core[first_key] = fene_bond_ia;
  BOOST_TEST(bond_map.at(first_key) == fene_bond_ia);
  BOOST_TEST(mock_core.at(first_key) == fene_bond_ia);
  BOOST_REQUIRE_EQUAL(bond_map.size(), 1);
  BOOST_REQUIRE_EQUAL(mock_core.size(), 1);
  // insert second element
  auto const acos_bond = AngleCosineBond(4.5, 3.14);
  auto const acos_bond_ia = std::make_shared<Bonded_IA_Parameters>(acos_bond);
  auto second_key = bond_map.insert(acos_bond_ia);
  BOOST_TEST(first_key != second_key);
  mock_core[second_key] = acos_bond_ia;
  BOOST_TEST(bond_map.at(second_key) == acos_bond_ia);
  BOOST_TEST(mock_core.at(second_key) == acos_bond_ia);
  BOOST_REQUIRE_EQUAL(bond_map.size(), 2);
  BOOST_REQUIRE_EQUAL(mock_core.size(), 2);
  // try to access non-existent element
  BOOST_CHECK_THROW(bond_map.at(bond_map.get_next_key()), std::out_of_range);
  BOOST_CHECK_THROW(mock_core.at(bond_map.get_next_key()), std::out_of_range);
  // check that both elements are "contained" in the map
  BOOST_TEST(bond_map.contains(first_key));
  BOOST_TEST(bond_map.contains(second_key));

  BOOST_REQUIRE_EQUAL(bond_map.count(first_key), 1);
  BOOST_REQUIRE_EQUAL(mock_core.count(first_key), 1);
  BOOST_REQUIRE_EQUAL(bond_map.count(second_key), 1);
  BOOST_REQUIRE_EQUAL(mock_core.count(second_key), 1);
  BOOST_REQUIRE_EQUAL(bond_map.count(bond_map.get_next_key()), 0);
  BOOST_REQUIRE_EQUAL(mock_core.count(bond_map.get_next_key()), 0);

  // delete an element
  bond_map.erase(first_key);
  mock_core.erase(first_key);
  BOOST_TEST(!bond_map.contains(first_key));
  BOOST_REQUIRE_EQUAL(bond_map.size(), 1);
  BOOST_REQUIRE_EQUAL(mock_core.size(), 1);
}
