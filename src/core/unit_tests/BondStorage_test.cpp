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

#define BOOST_TEST_MODULE BondStorage
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "BondStorage.hpp"

#include <algorithm>
#include <cstddef>
#include <span>
#include <utility>
#include <vector>

namespace {
using Tuple = std::pair<int, std::vector<int>>;

auto as_tuple(BondStorage::Bond const &bond) {
  return Tuple{bond.bond_id,
               {bond.particle_ids.begin(), bond.particle_ids.end()}};
}

/** @brief All bonds, sorted. */
auto all_bonds(BondStorage const &storage) {
  std::vector<Tuple> bonds;
  storage.for_each_bond(
      [&](BondStorage::Bond const &bond) { bonds.push_back(as_tuple(bond)); });
  std::ranges::sort(bonds);
  return bonds;
}

/** @brief Bonds of a particle, in the order they are visited. */
auto bonds_of(BondStorage const &storage, int pid) {
  std::vector<Tuple> bonds;
  storage.for_each_bond_of(pid, [&](BondStorage::Bond const &bond) {
    bonds.push_back(as_tuple(bond));
  });
  return bonds;
}

void add(BondStorage &storage, int bond_id, std::vector<int> const &ids) {
  storage.add(bond_id, ids);
}

bool remove(BondStorage &storage, int bond_id, std::vector<int> const &ids) {
  return storage.remove(bond_id, ids);
}
} // namespace

BOOST_AUTO_TEST_CASE(add_and_query) {
  BondStorage storage;
  BOOST_CHECK(storage.empty());

  add(storage, 0, {1, 2});
  add(storage, 1, {2, 1, 3});
  add(storage, 2, {4, 3, 2, 1});
  add(storage, 3, {5}); // bond without partners
  BOOST_CHECK_EQUAL(storage.size(), 4u);

  /* every bond is stored once, with the participants in creation order */
  BOOST_CHECK((all_bonds(storage) ==
               std::vector<Tuple>{
                   {0, {1, 2}}, {1, {2, 1, 3}}, {2, {4, 3, 2, 1}}, {3, {5}}}));

  /* bonds are found from any participant, in insertion order */
  BOOST_CHECK(
      (bonds_of(storage, 1) ==
       std::vector<Tuple>{{0, {1, 2}}, {1, {2, 1, 3}}, {2, {4, 3, 2, 1}}}));
  BOOST_CHECK((bonds_of(storage, 3) ==
               std::vector<Tuple>{{1, {2, 1, 3}}, {2, {4, 3, 2, 1}}}));
  BOOST_CHECK((bonds_of(storage, 4) == std::vector<Tuple>{{2, {4, 3, 2, 1}}}));
  BOOST_CHECK((bonds_of(storage, 5) == std::vector<Tuple>{{3, {5}}}));
  BOOST_CHECK(bonds_of(storage, 6).empty());

  /* accessors */
  storage.for_each_bond_of(3, [](BondStorage::Bond const &bond) {
    BOOST_CHECK(bond.has_participant(3));
    BOOST_CHECK(not bond.has_participant(5));
    BOOST_CHECK_EQUAL(bond.owner_id(), bond.particle_ids.front());
    BOOST_CHECK_EQUAL(bond.partner_ids().size(), bond.particle_ids.size() - 1u);
    BOOST_CHECK_EQUAL(bond.partner_ids().front(), bond.particle_ids[1]);
  });

  storage.clear();
  BOOST_CHECK(storage.empty());
  BOOST_CHECK(all_bonds(storage).empty());
  BOOST_CHECK(bonds_of(storage, 1).empty());
}

BOOST_AUTO_TEST_CASE(remove_matching) {
  BondStorage storage;
  add(storage, 0, {1, 2});
  add(storage, 1, {2, 1, 3});

  /* no match: wrong bond id, wrong participants, wrong arity, no ids */
  BOOST_CHECK(not remove(storage, 1, {1, 2}));
  BOOST_CHECK(not remove(storage, 0, {1, 3}));
  BOOST_CHECK(not remove(storage, 0, {7, 8}));
  BOOST_CHECK(not remove(storage, 1, {2, 1}));
  BOOST_CHECK(not remove(storage, 0, {}));
  BOOST_CHECK(not remove(storage, 0, {1, 2, 3, 4, 5}));
  BOOST_CHECK_EQUAL(storage.size(), 2u);

  /* a bond can be removed from any participant, in any order */
  BOOST_CHECK(remove(storage, 1, {3, 2, 1}));
  BOOST_CHECK(remove(storage, 0, {2, 1}));
  BOOST_CHECK(storage.empty());
  BOOST_CHECK(bonds_of(storage, 1).empty());
  BOOST_CHECK(bonds_of(storage, 2).empty());
  BOOST_CHECK(bonds_of(storage, 3).empty());

  /* a removed bond cannot be removed again */
  BOOST_CHECK(not remove(storage, 0, {1, 2}));
}

BOOST_AUTO_TEST_CASE(remove_duplicates) {
  BondStorage storage;
  /* the same bond, stored twice with swapped ownership and once more */
  add(storage, 0, {1, 2});
  add(storage, 0, {2, 1});
  add(storage, 0, {1, 2});
  BOOST_CHECK_EQUAL(storage.size(), 3u);

  /* the bond owned by the first particle is preferred */
  BOOST_CHECK(remove(storage, 0, {2, 1}));
  BOOST_CHECK(
      (all_bonds(storage) == std::vector<Tuple>{{0, {1, 2}}, {0, {1, 2}}}));
  /* otherwise any bond with the same participants is removed */
  BOOST_CHECK(remove(storage, 0, {2, 1}));
  BOOST_CHECK((all_bonds(storage) == std::vector<Tuple>{{0, {1, 2}}}));
  BOOST_CHECK(remove(storage, 0, {1, 2}));
  BOOST_CHECK(storage.empty());
}

BOOST_AUTO_TEST_CASE(slot_reuse) {
  BondStorage storage;
  add(storage, 0, {1, 2});
  add(storage, 0, {2, 3});
  add(storage, 0, {3, 4});
  BOOST_CHECK(remove(storage, 0, {2, 3}));
  add(storage, 1, {5, 3});
  add(storage, 1, {6, 3});
  BOOST_CHECK_EQUAL(storage.size(), 4u);
  BOOST_CHECK(
      (all_bonds(storage) ==
       std::vector<Tuple>{{0, {1, 2}}, {0, {3, 4}}, {1, {5, 3}}, {1, {6, 3}}}));
  /* the bonds of a particle keep their insertion order */
  BOOST_CHECK((bonds_of(storage, 3) ==
               std::vector<Tuple>{{0, {3, 4}}, {1, {5, 3}}, {1, {6, 3}}}));
  BOOST_CHECK((bonds_of(storage, 2) == std::vector<Tuple>{{0, {1, 2}}}));
}

BOOST_AUTO_TEST_CASE(repeated_participant) {
  BondStorage storage;
  /* a particle listed twice in a bond is visited once */
  add(storage, 0, {1, 1});
  add(storage, 1, {2, 1, 1});
  BOOST_CHECK((bonds_of(storage, 1) ==
               std::vector<Tuple>{{0, {1, 1}}, {1, {2, 1, 1}}}));
  BOOST_CHECK(remove(storage, 1, {1, 1, 2}));
  BOOST_CHECK(remove(storage, 0, {1, 1}));
  BOOST_CHECK(storage.empty());
  BOOST_CHECK(bonds_of(storage, 1).empty());
}

BOOST_AUTO_TEST_CASE(owned_bonds) {
  /* bonds owned by a particle, in the order they are visited */
  auto const owned_by = [](BondStorage const &storage, int pid) {
    std::vector<Tuple> bonds;
    storage.owned_bonds().for_each_bond_of(
        pid, [&](int bond_id, std::span<const int> partner_ids) {
          bonds.push_back({bond_id, {partner_ids.begin(), partner_ids.end()}});
        });
    return bonds;
  };
  auto const counts = [](BondStorage const &storage, int pid) {
    auto const c = storage.owned_bonds().counts(pid);
    return std::vector<int>{c.pair, c.angle, c.dihedral};
  };

  BondStorage storage;
  /* empty storage, unknown and invalid particle ids */
  BOOST_CHECK(owned_by(storage, 0).empty());
  BOOST_CHECK(owned_by(storage, -1).empty());
  BOOST_CHECK((counts(storage, 0) == std::vector<int>{0, 0, 0}));

  add(storage, 0, {2, 1});
  add(storage, 1, {2, 1, 3});
  add(storage, 2, {4, 3, 2, 1});
  add(storage, 3, {5}); // bonds without partners are not included
  add(storage, 0, {2, 7});
  add(storage, 4, {2, 6, 5});
  add(storage, 0, {1, 2});

  /* only the bonds a particle owns, grouped by number of partners, with
   * the partners in creation order */
  BOOST_CHECK(
      (owned_by(storage, 2) ==
       std::vector<Tuple>{{0, {1}}, {0, {7}}, {1, {1, 3}}, {4, {6, 5}}}));
  BOOST_CHECK((owned_by(storage, 1) == std::vector<Tuple>{{0, {2}}}));
  BOOST_CHECK((owned_by(storage, 4) == std::vector<Tuple>{{2, {3, 2, 1}}}));
  BOOST_CHECK(owned_by(storage, 3).empty());
  BOOST_CHECK(owned_by(storage, 5).empty());
  BOOST_CHECK(owned_by(storage, 7).empty());
  BOOST_CHECK(owned_by(storage, 100).empty());
  BOOST_CHECK((counts(storage, 2) == std::vector<int>{2, 2, 0}));
  BOOST_CHECK((counts(storage, 4) == std::vector<int>{0, 0, 1}));
  BOOST_CHECK((counts(storage, 5) == std::vector<int>{0, 0, 0}));
  BOOST_CHECK((counts(storage, 100) == std::vector<int>{0, 0, 0}));

  /* over all particles, every bond with partners is visited once */
  std::size_t n_bonds = 0u;
  for (int pid = 0; pid < 10; ++pid) {
    n_bonds += owned_by(storage, pid).size();
  }
  BOOST_CHECK_EQUAL(n_bonds, storage.size() - 1u);

  /* the view follows modifications of the storage */
  BOOST_CHECK(remove(storage, 0, {2, 1}));
  BOOST_CHECK((owned_by(storage, 2) ==
               std::vector<Tuple>{{0, {7}}, {1, {1, 3}}, {4, {6, 5}}}));
  BOOST_CHECK((owned_by(storage, 1) == std::vector<Tuple>{{0, {2}}}));
  add(storage, 5, {9, 2});
  BOOST_CHECK((owned_by(storage, 9) == std::vector<Tuple>{{5, {2}}}));
  storage.clear();
  BOOST_CHECK(owned_by(storage, 2).empty());
  BOOST_CHECK((counts(storage, 2) == std::vector<int>{0, 0, 0}));
}
