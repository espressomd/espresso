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

#define BOOST_TEST_MODULE BondStore
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "BondStore.hpp"

#include <cstddef>
#include <span>
#include <utility>
#include <vector>

namespace {
using Tuple = std::pair<int, std::vector<int>>;

auto as_tuples(BondRecordRange const &range) {
  std::vector<Tuple> out;
  for (auto const bond : range) {
    auto const ids = bond.particle_ids();
    out.push_back({bond.bond_id(), {ids.begin(), ids.end()}});
  }
  return out;
}

auto counts(BondStore const &store, int pid) {
  auto const c = store.owned_counts(pid);
  return std::vector<int>{c.pair, c.angle, c.dihedral};
}

void add(BondStore &store, int pid, int bond_id, std::vector<int> const &ids) {
  store.add(pid, bond_id, ids);
}

bool remove(BondStore &store, int pid, int bond_id,
            std::vector<int> const &ids) {
  return store.remove(pid, bond_id, ids);
}
} // namespace

BOOST_AUTO_TEST_CASE(add_and_view) {
  BondStore store;
  BOOST_CHECK(store.empty());
  BOOST_CHECK(store.view(2).empty());
  BOOST_CHECK((counts(store, 2) == std::vector<int>{0, 0, 0}));

  /* records of particle 2: owned and not owned, any arity */
  add(store, 2, 0, {2, 1});
  add(store, 2, 1, {1, 2, 3});
  add(store, 2, 2, {4, 3, 2, 1});
  add(store, 2, 3, {2}); // bond without partners: owned, but not counted
  add(store, 2, 4, {2, 6, 5});
  /* the owned records come first, both groups in insertion order */
  BOOST_CHECK(
      (as_tuples(store.view(2)) == std::vector<Tuple>{{0, {2, 1}},
                                                      {3, {2}},
                                                      {4, {2, 6, 5}},
                                                      {1, {1, 2, 3}},
                                                      {2, {4, 3, 2, 1}}}));
  BOOST_CHECK((as_tuples(store.owned(2)) ==
               std::vector<Tuple>{{0, {2, 1}}, {3, {2}}, {4, {2, 6, 5}}}));
  BOOST_CHECK((counts(store, 2) == std::vector<int>{1, 1, 0}));
  BOOST_CHECK_EQUAL(store.view(2).size(), 5u);
  BOOST_CHECK(store.contains(2, 1, std::vector<int>{1, 2, 3}));
  BOOST_CHECK(not store.contains(2, 1, std::vector<int>{3, 2, 1}));
  BOOST_CHECK(not store.contains(7, 1, std::vector<int>{1, 2, 3}));

  /* record accessors */
  auto const bond = *store.view(2).begin();
  BOOST_CHECK_EQUAL(bond.bond_id(), 0);
  BOOST_CHECK_EQUAL(bond.owner_id(), 2);
  BOOST_CHECK_EQUAL(bond.partner_ids().size(), 1u);
  BOOST_CHECK_EQUAL(bond.partner_ids().front(), 1);
  BOOST_CHECK(bond.has_participant(1));
  BOOST_CHECK(not bond.has_participant(3));
  BOOST_CHECK(bond.matches(0, std::vector<int>{2, 1}));
  BOOST_CHECK(not bond.matches(0, std::vector<int>{1, 2}));

  /* other particles are independent */
  add(store, 1, 0, {2, 1});
  BOOST_CHECK_EQUAL(store.n_particles(), 2u);
  BOOST_CHECK(store.owned(1).empty());
  BOOST_CHECK_EQUAL(store.view(1).size(), 1u);
  BOOST_CHECK_EQUAL(store.view(2).size(), 5u);

  store.erase(2);
  BOOST_CHECK(store.view(2).empty());
  BOOST_CHECK_EQUAL(store.view(1).size(), 1u);
  store.clear();
  BOOST_CHECK(store.empty());
  BOOST_CHECK(store.view(1).empty());
}

BOOST_AUTO_TEST_CASE(remove_exact_order) {
  BondStore store;
  /* two dihedrals and two angles with the same participants, but
   * different partner orders, plus a bond stored twice */
  add(store, 1, 0, {1, 2, 3, 4});
  add(store, 1, 0, {1, 3, 2, 4});
  add(store, 1, 1, {1, 2, 3});
  add(store, 1, 1, {1, 3, 2});
  add(store, 1, 2, {1, 2});
  add(store, 1, 2, {1, 2});
  BOOST_CHECK((counts(store, 1) == std::vector<int>{2, 2, 2}));

  /* no match: wrong bond id, participants, order, arity, no ids, particle */
  BOOST_CHECK(not remove(store, 1, 1, {1, 2}));
  BOOST_CHECK(not remove(store, 1, 2, {2, 1}));
  BOOST_CHECK(not remove(store, 1, 0, {1, 4, 3, 2}));
  BOOST_CHECK(not remove(store, 1, 0, {1, 2, 3}));
  BOOST_CHECK(not remove(store, 1, 0, {}));
  BOOST_CHECK(not remove(store, 7, 2, {1, 2}));
  BOOST_CHECK_EQUAL(store.view(1).size(), 6u);

  /* the record with the given order is removed, not the first stored */
  BOOST_CHECK(remove(store, 1, 0, {1, 3, 2, 4}));
  BOOST_CHECK(remove(store, 1, 1, {1, 3, 2}));
  /* a bond stored twice is removed one occurrence at a time */
  BOOST_CHECK(remove(store, 1, 2, {1, 2}));
  BOOST_CHECK(
      (as_tuples(store.view(1)) ==
       std::vector<Tuple>{{0, {1, 2, 3, 4}}, {1, {1, 2, 3}}, {2, {1, 2}}}));
  BOOST_CHECK((counts(store, 1) == std::vector<int>{1, 1, 1}));
  BOOST_CHECK(remove(store, 1, 2, {1, 2}));
  BOOST_CHECK(not remove(store, 1, 2, {1, 2}));
  BOOST_CHECK(remove(store, 1, 0, {1, 2, 3, 4}));
  BOOST_CHECK(remove(store, 1, 1, {1, 2, 3}));
  /* a particle without records is forgotten */
  BOOST_CHECK(store.empty());
}

BOOST_AUTO_TEST_CASE(owned_prefix_follows_changes) {
  BondStore store;
  add(store, 1, 0, {1, 2});
  add(store, 1, 1, {3, 1, 4});
  add(store, 1, 0, {5, 1});
  add(store, 1, 2, {1, 6, 7, 8});
  add(store, 1, 0, {1, 9});
  BOOST_CHECK(
      (as_tuples(store.owned(1)) ==
       std::vector<Tuple>{{0, {1, 2}}, {2, {1, 6, 7, 8}}, {0, {1, 9}}}));

  BOOST_CHECK(remove(store, 1, 2, {1, 6, 7, 8}));
  BOOST_CHECK(remove(store, 1, 1, {3, 1, 4}));
  add(store, 1, 3, {1, 10});
  add(store, 1, 3, {11, 1});
  BOOST_CHECK(
      (as_tuples(store.view(1)) ==
       std::vector<Tuple>{
           {0, {1, 2}}, {0, {1, 9}}, {3, {1, 10}}, {0, {5, 1}}, {3, {11, 1}}}));
  BOOST_CHECK((counts(store, 1) == std::vector<int>{3, 0, 0}));

  /* remove all records with particle 5 or 9 */
  auto const n = store.remove_if(1, [](BondRecord const &bond) {
    return bond.has_participant(5) or bond.has_participant(9);
  });
  BOOST_CHECK_EQUAL(n, 2u);
  BOOST_CHECK((as_tuples(store.view(1)) ==
               std::vector<Tuple>{{0, {1, 2}}, {3, {1, 10}}, {3, {11, 1}}}));
  BOOST_CHECK((as_tuples(store.owned(1)) ==
               std::vector<Tuple>{{0, {1, 2}}, {3, {1, 10}}}));
  BOOST_CHECK((counts(store, 1) == std::vector<int>{2, 0, 0}));
  BOOST_CHECK_EQUAL(
      store.remove_if(1, [](BondRecord const &) { return false; }), 0u);
  BOOST_CHECK_EQUAL(store.remove_if(1, [](BondRecord const &) { return true; }),
                    3u);
  BOOST_CHECK(store.empty());
}

BOOST_AUTO_TEST_CASE(compaction_and_index_access) {
  BondStore store;
  /* interleaved changes relocate the blocks */
  for (int i = 0; i < 50; ++i) {
    for (int pid = 0; pid < 4; ++pid) {
      add(store, pid, i, {pid, 100 + i});
    }
  }
  BOOST_CHECK(remove(store, 2, 7, {2, 107}));
  add(store, 3, 99, {5, 3}); // not owned
  BOOST_CHECK(not store.is_indexed());
  auto const before = std::vector<std::vector<Tuple>>{
      as_tuples(store.view(0)), as_tuples(store.view(1)),
      as_tuples(store.view(2)), as_tuples(store.view(3))};

  /* compaction in a given order, with a particle without records and a
   * particle with records that is not in the order */
  auto const order = std::vector<int>{3, 42, 1, 2};
  store.compact(order);
  BOOST_CHECK(store.is_indexed());
  /* the buffer holds exactly the records */
  auto n_ints = std::size_t{0u};
  for (int pid = 0; pid < 4; ++pid) {
    n_ints += store.view(pid).data().size();
    BOOST_CHECK((as_tuples(store.view(pid)) == before[pid]));
  }
  BOOST_CHECK_EQUAL(store.capacity(), n_ints);
  BOOST_CHECK_EQUAL(store.size_in_use(), n_ints);
  BOOST_CHECK(not store.is_wasteful());
  /* the blocks are contiguous, in the given order */
  BOOST_CHECK(store.view(1).data().data() ==
              store.view(3).data().data() + store.view(3).data().size());
  BOOST_CHECK(store.view(2).data().data() ==
              store.view(1).data().data() + store.view(1).data().size());
  /* access by position in the order */
  BOOST_CHECK((as_tuples(store.owned_at(0)) == as_tuples(store.owned(3))));
  BOOST_CHECK(store.owned_at(1).empty());
  BOOST_CHECK((as_tuples(store.owned_at(2)) == as_tuples(store.owned(1))));
  BOOST_CHECK_EQUAL(store.owned_counts_at(3).pair, 49);
  BOOST_CHECK_EQUAL(store.owned_counts_at(0).pair, 50);
  /* the access by position can also be set up without moving records */
  store.index_resize(3u);
  store.set_index(0u, 2);
  store.set_index(1u, 42);
  store.set_index(2u, 0);
  store.index_done();
  BOOST_CHECK(store.is_indexed());
  BOOST_CHECK((as_tuples(store.owned_at(0)) == as_tuples(store.owned(2))));
  BOOST_CHECK(store.owned_at(1).empty());
  BOOST_CHECK_EQUAL(store.owned_counts_at(2).pair, 50);
  /* any change invalidates the access by position */
  add(store, 1, 0, {1, 2});
  BOOST_CHECK(not store.is_indexed());
}

BOOST_AUTO_TEST_CASE(extract_and_insert) {
  BondStore source;
  add(source, 2, 0, {2, 1});
  add(source, 2, 1, {1, 2, 3});
  add(source, 2, 2, {4, 3, 2, 1});
  add(source, 5, 0, {5, 6});
  auto const records_2 = as_tuples(source.view(2));
  auto const records_5 = as_tuples(source.view(5));

  /* particles leave in the order 5, 7 (no records), 2 */
  std::vector<int> buffer;
  source.extract(5, buffer);
  source.extract(7, buffer);
  source.extract(2, buffer);
  BOOST_CHECK(source.empty());

  BondStore target;
  add(target, 2, 9, {9, 2}); // replaced
  auto rest = std::span<const int>(buffer);
  rest = target.insert(5, rest);
  rest = target.insert(7, rest);
  rest = target.insert(2, rest);
  BOOST_CHECK(rest.empty());
  BOOST_CHECK((as_tuples(target.view(2)) == records_2));
  BOOST_CHECK((as_tuples(target.view(5)) == records_5));
  BOOST_CHECK(target.view(7).empty());
  BOOST_CHECK((counts(target, 2) == std::vector<int>{1, 0, 0}));
  BOOST_CHECK((as_tuples(target.owned(2)) == std::vector<Tuple>{{0, {2, 1}}}));

  /* a copy does not remove the records */
  std::vector<int> copy;
  target.copy_to(2, copy);
  BOOST_CHECK_EQUAL(target.view(2).size(), 3u);
  BondStore other;
  other.insert(2, copy);
  BOOST_CHECK((as_tuples(other.view(2)) == records_2));
}

BOOST_AUTO_TEST_CASE(many_bonds) {
  /* the center of a star polymer */
  BondStore store;
  int const n_arms = 40000;
  for (int i = 1; i <= n_arms; ++i) {
    store.add(0, 0, std::vector<int>{0, i});
    store.add(i, 0, std::vector<int>{0, i});
  }
  BOOST_CHECK_EQUAL(store.view(0).size(), static_cast<std::size_t>(n_arms));
  BOOST_CHECK_EQUAL(store.owned_counts(0).pair, n_arms);
  BOOST_CHECK(store.contains(0, 0, std::vector<int>{0, n_arms}));
  BOOST_CHECK(store.remove(0, 0, std::vector<int>{0, 1}));
  BOOST_CHECK_EQUAL(store.owned_counts(0).pair, n_arms - 1);
  BOOST_CHECK_EQUAL(store.view(17).size(), 1u);
  /* the unused ints are reclaimed while the hub keeps growing */
  BOOST_CHECK_LT(store.capacity(), 4u * (2u * 4u * n_arms));
}
