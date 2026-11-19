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
  storage.for_each_bond_with(pid, [&](BondStorage::Bond const &bond) {
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

  /* bonds are found from any participant */
  BOOST_CHECK(
      (bonds_of(storage, 1) ==
       std::vector<Tuple>{{0, {1, 2}}, {1, {2, 1, 3}}, {2, {4, 3, 2, 1}}}));
  BOOST_CHECK((bonds_of(storage, 3) ==
               std::vector<Tuple>{{1, {2, 1, 3}}, {2, {4, 3, 2, 1}}}));
  BOOST_CHECK((bonds_of(storage, 4) == std::vector<Tuple>{{2, {4, 3, 2, 1}}}));
  BOOST_CHECK((bonds_of(storage, 5) == std::vector<Tuple>{{3, {5}}}));
  BOOST_CHECK(bonds_of(storage, 6).empty());

  /* accessors */
  storage.for_each_bond_with(3, [](BondStorage::Bond const &bond) {
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

  /* a bond can be removed in exact order */
  BOOST_CHECK(not remove(storage, 1, {3, 2, 1}));
  BOOST_CHECK(not remove(storage, 0, {2, 1}));
  BOOST_CHECK(remove(storage, 1, {2, 1, 3}));
  BOOST_CHECK(remove(storage, 0, {1, 2}));
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

  /* only the exact tuple is removed: (2, 1) is a different bond than (1, 2) */
  BOOST_CHECK(remove(storage, 0, {2, 1}));
  BOOST_CHECK(
      (all_bonds(storage) == std::vector<Tuple>{{0, {1, 2}}, {0, {1, 2}}}));
  BOOST_CHECK(not remove(storage, 0, {2, 1}));
  /* a bond stored twice is removed one occurrence at a time */
  BOOST_CHECK(remove(storage, 0, {1, 2}));
  BOOST_CHECK((all_bonds(storage) == std::vector<Tuple>{{0, {1, 2}}}));
  BOOST_CHECK(remove(storage, 0, {1, 2}));
  BOOST_CHECK(storage.empty());
}

BOOST_AUTO_TEST_CASE(remove_exact_order) {
  BondStorage storage;
  /* two dihedrals with the same participants, but different bonds */
  add(storage, 0, {1, 2, 3, 4});
  add(storage, 0, {1, 3, 2, 4});
  /* two angles with the same participants, but different vertices */
  add(storage, 1, {1, 2, 3});
  add(storage, 1, {1, 3, 2});

  /* the bond with the given partner order is removed, not the one stored
   * first */
  BOOST_CHECK(remove(storage, 0, {1, 3, 2, 4}));
  BOOST_CHECK(remove(storage, 1, {1, 3, 2}));
  BOOST_CHECK((all_bonds(storage) ==
               std::vector<Tuple>{{0, {1, 2, 3, 4}}, {1, {1, 2, 3}}}));
  BOOST_CHECK((storage.owned_counts(1).angle == 1));
  BOOST_CHECK((storage.owned_counts(1).dihedral == 1));

  /* a different partner order is a different bond: nothing is removed */
  BOOST_CHECK(not remove(storage, 0, {1, 4, 3, 2}));
  BOOST_CHECK(not remove(storage, 1, {3, 2, 1}));
  BOOST_CHECK_EQUAL(storage.size(), 2u);
  BOOST_CHECK(remove(storage, 0, {1, 2, 3, 4}));
  BOOST_CHECK(remove(storage, 1, {1, 2, 3}));
  BOOST_CHECK(storage.empty());
}

BOOST_AUTO_TEST_CASE(visiting_order) {
  BondStorage storage;
  add(storage, 0, {2, 1});
  add(storage, 1, {1, 2});
  add(storage, 2, {3, 1, 2});
  add(storage, 3, {2, 1});
  add(storage, 4, {1, 3});
  /* the bonds a particle owns come first; both groups are in insertion
   * order, the bonds of other owners being grouped by owner */
  BOOST_CHECK((
      bonds_of(storage, 1) ==
      std::vector<Tuple>{
          {1, {1, 2}}, {4, {1, 3}}, {0, {2, 1}}, {3, {2, 1}}, {2, {3, 1, 2}}}));
  /* a particle is found from the remaining bonds of an owner */
  BOOST_CHECK(remove(storage, 0, {2, 1}));
  BOOST_CHECK((bonds_of(storage, 1) ==
               std::vector<Tuple>{
                   {1, {1, 2}}, {4, {1, 3}}, {3, {2, 1}}, {2, {3, 1, 2}}}));
  BOOST_CHECK(remove(storage, 3, {2, 1}));
  BOOST_CHECK((bonds_of(storage, 1) ==
               std::vector<Tuple>{{1, {1, 2}}, {4, {1, 3}}, {2, {3, 1, 2}}}));
  /* negative ids are never found */
  BOOST_CHECK(bonds_of(storage, -1).empty());
  BOOST_CHECK(not remove(storage, 1, {-1, 2}));
}

BOOST_AUTO_TEST_CASE(removal_and_reuse) {
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
  /* the remaining bonds of a particle keep their order */
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
  BOOST_CHECK(not remove(storage, 1, {1, 1, 2}));
  BOOST_CHECK(remove(storage, 1, {2, 1, 1}));
  BOOST_CHECK(remove(storage, 0, {1, 1}));
  BOOST_CHECK(storage.empty());
  BOOST_CHECK(bonds_of(storage, 1).empty());
}

BOOST_AUTO_TEST_CASE(owned_bonds) {
  /* bonds owned by a particle, in the order they are visited */
  auto const owned_by = [](BondStorage const &storage, int pid) {
    std::vector<Tuple> bonds;
    storage.for_each_owned_bond(pid, [&](BondStorage::Bond const &bond) {
      auto const partner_ids = bond.partner_ids();
      bonds.push_back({bond.bond_id, {partner_ids.begin(), partner_ids.end()}});
    });
    return bonds;
  };
  auto const counts = [](BondStorage const &storage, int pid) {
    auto const c = storage.owned_counts(pid);
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
  add(storage, 3, {5}); // bond without partners: owned, but not counted
  add(storage, 0, {2, 7});
  add(storage, 4, {2, 6, 5});
  add(storage, 0, {1, 2});

  /* only the bonds a particle owns, in insertion order, with the partners
   * in creation order */
  BOOST_CHECK(
      (owned_by(storage, 2) ==
       std::vector<Tuple>{{0, {1}}, {1, {1, 3}}, {0, {7}}, {4, {6, 5}}}));
  BOOST_CHECK((owned_by(storage, 1) == std::vector<Tuple>{{0, {2}}}));
  BOOST_CHECK((owned_by(storage, 4) == std::vector<Tuple>{{2, {3, 2, 1}}}));
  BOOST_CHECK(owned_by(storage, 3).empty());
  BOOST_CHECK((owned_by(storage, 5) == std::vector<Tuple>{{3, {}}}));
  BOOST_CHECK(owned_by(storage, 7).empty());
  BOOST_CHECK(owned_by(storage, 100).empty());
  BOOST_CHECK((counts(storage, 2) == std::vector<int>{2, 2, 0}));
  BOOST_CHECK((counts(storage, 4) == std::vector<int>{0, 0, 1}));
  BOOST_CHECK((counts(storage, 5) == std::vector<int>{0, 0, 0}));
  BOOST_CHECK((counts(storage, 100) == std::vector<int>{0, 0, 0}));

  /* over all particles, every bond is visited once */
  std::size_t n_bonds = 0u;
  for (int pid = 0; pid < 10; ++pid) {
    n_bonds += owned_by(storage, pid).size();
  }
  BOOST_CHECK_EQUAL(n_bonds, storage.size());

  /* any_owned_bond() only considers the bonds a particle owns */
  auto const has_id = [](int bond_id) {
    return
        [=](BondStorage::Bond const &bond) { return bond.bond_id == bond_id; };
  };
  BOOST_CHECK(storage.any_owned_bond(2, has_id(1)));
  BOOST_CHECK(not storage.any_owned_bond(1, has_id(1)));
  BOOST_CHECK(storage.any_bond_with(1, has_id(1)));
  BOOST_CHECK(not storage.any_owned_bond(100, has_id(1)));
  BOOST_CHECK(not storage.any_owned_bond(-1, has_id(1)));

  /* the view follows modifications of the storage */
  BOOST_CHECK(remove(storage, 0, {2, 1}));
  BOOST_CHECK((owned_by(storage, 2) ==
               std::vector<Tuple>{{1, {1, 3}}, {0, {7}}, {4, {6, 5}}}));
  BOOST_CHECK((counts(storage, 2) == std::vector<int>{1, 2, 0}));
  BOOST_CHECK((owned_by(storage, 1) == std::vector<Tuple>{{0, {2}}}));
  add(storage, 5, {9, 2});
  BOOST_CHECK((owned_by(storage, 9) == std::vector<Tuple>{{5, {2}}}));
  storage.clear();
  BOOST_CHECK(owned_by(storage, 2).empty());
  BOOST_CHECK((counts(storage, 2) == std::vector<int>{0, 0, 0}));
}
