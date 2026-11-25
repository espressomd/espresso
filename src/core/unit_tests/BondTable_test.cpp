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

#define BOOST_TEST_MODULE BondTable
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "BondTable.hpp"

#include <cstddef>
#include <iterator>
#include <ranges>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

namespace {
using Tuple = std::pair<int, std::vector<int>>;

template <typename Range> auto as_tuples(Range &&range) {
  std::vector<Tuple> out;
  for (auto const bond : range) {
    auto const ids = bond.particle_ids();
    out.push_back({bond.bond_id(), {ids.begin(), ids.end()}});
  }
  return out;
}

/** @brief Add a bond whose participants in @p local reside here. */
void add(BondTable &table, int bond_id, std::vector<int> const &ids,
         std::vector<int> const &local) {
  table.add(bond_id, ids, local);
}

bool remove(BondTable &table, int bond_id, std::vector<int> const &ids) {
  return table.remove(bond_id, ids);
}

/** @brief Number of non-free rows. */
auto n_live_rows(BondTable const &table) {
  return static_cast<std::size_t>(std::ranges::count_if(
      table.rows(), [](BondRow const &row) { return not row.is_free(); }));
}
} // namespace

BOOST_AUTO_TEST_CASE(one_row_per_bond) {
  BondTable table;
  BOOST_CHECK(table.empty());
  BOOST_CHECK(table.view(1).empty());

  /* particles 1, 2, 3 reside here, 4 does not */
  add(table, 0, {1, 2}, {1, 2});
  add(table, 1, {2, 1, 3}, {2, 1, 3});
  add(table, 2, {4, 3, 2, 1}, {3, 2, 1});
  add(table, 3, {2}, {2}); // bond without partners
  /* every bond is stored once, whatever its number of local participants */
  BOOST_CHECK_EQUAL(table.size(), 4u);
  BOOST_CHECK_EQUAL(n_live_rows(table), 4u);
  BOOST_CHECK_EQUAL(table.n_particles(), 3u);

  /* the bonds are found from every local participant, in insertion order */
  BOOST_CHECK((as_tuples(table.view(2)) ==
               std::vector<Tuple>{
                   {0, {1, 2}}, {1, {2, 1, 3}}, {2, {4, 3, 2, 1}}, {3, {2}}}));
  BOOST_CHECK((as_tuples(table.view(3)) ==
               std::vector<Tuple>{{1, {2, 1, 3}}, {2, {4, 3, 2, 1}}}));
  BOOST_CHECK(table.view(4).empty());
  /* owned bonds only */
  BOOST_CHECK((as_tuples(table.owned(2)) ==
               std::vector<Tuple>{{1, {2, 1, 3}}, {3, {2}}}));
  BOOST_CHECK(table.owned(3).empty());
  BOOST_CHECK(table.owned(4).empty());

  /* a bond is identified by its id and its participants in order */
  BOOST_CHECK(table.contains(1, std::vector<int>{2, 1, 3}));
  BOOST_CHECK(not table.contains(1, std::vector<int>{2, 3, 1}));
  BOOST_CHECK(not table.contains(0, std::vector<int>{2, 1}));
  /* a stored bond is not stored again, a reversed one is a different bond */
  add(table, 0, {1, 2}, {1, 2});
  BOOST_CHECK_EQUAL(table.size(), 4u);
  add(table, 0, {2, 1}, {2, 1});
  BOOST_CHECK_EQUAL(table.size(), 5u);
  /* a bond without local participant is not stored */
  add(table, 0, {7, 8}, {});
  BOOST_CHECK_EQUAL(table.size(), 5u);

  /* rows give access to the bonds for the force kernels */
  auto const &row = table.rows()[0];
  BOOST_CHECK_EQUAL(row.bond_id, 0);
  BOOST_CHECK_EQUAL(row.n, 2u);
  BOOST_CHECK_EQUAL(row.local_mask, 0b11u);
  BOOST_CHECK(row.is_local(1) and row.is_local(2));
  BOOST_CHECK_EQUAL(row.record().owner_id(), 1);
  BOOST_CHECK_EQUAL(row.record().partner_ids()[0], 2);

  table.clear();
  BOOST_CHECK(table.empty());
  BOOST_CHECK(table.view(2).empty());
}

BOOST_AUTO_TEST_CASE(remove_and_reuse) {
  BondTable table;
  add(table, 0, {1, 2, 3, 4}, {1, 2, 3, 4});
  add(table, 0, {1, 3, 2, 4}, {1, 3, 2, 4});
  add(table, 1, {1, 2}, {1, 2});

  /* no match: wrong bond id, participants, order, no ids */
  BOOST_CHECK(not remove(table, 1, {2, 1}));
  BOOST_CHECK(not remove(table, 0, {1, 4, 3, 2}));
  BOOST_CHECK(not remove(table, 0, {}));
  BOOST_CHECK(not remove(table, 1, {7, 8}));
  BOOST_CHECK_EQUAL(table.size(), 3u);

  /* the bond with the given order is removed, from all participants */
  BOOST_CHECK(remove(table, 0, {1, 3, 2, 4}));
  BOOST_CHECK(
      (as_tuples(table.view(3)) == std::vector<Tuple>{{0, {1, 2, 3, 4}}}));
  BOOST_CHECK(not remove(table, 0, {1, 3, 2, 4}));
  BOOST_CHECK_EQUAL(table.size(), 2u);
  /* the freed row is reused */
  auto const n_rows = table.rows().size();
  add(table, 2, {5, 6}, {5, 6});
  BOOST_CHECK_EQUAL(table.rows().size(), n_rows);
  BOOST_CHECK(remove(table, 0, {1, 2, 3, 4}));
  BOOST_CHECK(remove(table, 1, {1, 2}));
  BOOST_CHECK(remove(table, 2, {5, 6}));
  BOOST_CHECK(table.empty());
  BOOST_CHECK_EQUAL(table.n_particles(), 0u);
}

BOOST_AUTO_TEST_CASE(remove_if_and_erase) {
  BondTable table;
  add(table, 0, {1, 2}, {1, 2});
  add(table, 1, {3, 1, 4}, {3, 1});
  add(table, 0, {5, 1}, {1});
  add(table, 2, {1, 6, 7, 8}, {1});

  /* remove the bonds of particle 1 with particle 5 or bond id 1 */
  auto const n = table.remove_if(1, [](BondRecord const &bond) {
    return bond.has_participant(5) or bond.bond_id() == 1;
  });
  BOOST_CHECK_EQUAL(n, 2u);
  BOOST_CHECK((as_tuples(table.view(1)) ==
               std::vector<Tuple>{{0, {1, 2}}, {2, {1, 6, 7, 8}}}));
  /* the removed bonds are gone for the other participants too */
  BOOST_CHECK(table.view(3).empty());

  /* a particle that leaves no longer lists its bonds; bonds with another
   * local participant are kept, the others are freed */
  table.erase(1);
  BOOST_CHECK(table.view(1).empty());
  BOOST_CHECK((as_tuples(table.view(2)) == std::vector<Tuple>{{0, {1, 2}}}));
  BOOST_CHECK_EQUAL(table.size(), 1u);
}

BOOST_AUTO_TEST_CASE(remove_bonds_of) {
  /* particles 1, 2, 3 reside here, 4 and 5 do not */
  BondTable table;
  add(table, 0, {1, 2}, {1, 2});
  add(table, 0, {4, 1}, {1});
  add(table, 1, {2, 4, 3}, {2, 3});
  add(table, 0, {3, 5}, {3});
  BOOST_CHECK(not table.rows()[1].is_local(4));

  /* the bonds of a particle residing elsewhere are found too */
  BOOST_CHECK_EQUAL(table.remove_bonds_of(4), 2u);
  BOOST_CHECK((as_tuples(table.view(1)) == std::vector<Tuple>{{0, {1, 2}}}));
  BOOST_CHECK((as_tuples(table.view(3)) == std::vector<Tuple>{{0, {3, 5}}}));
  BOOST_CHECK_EQUAL(table.remove_bonds_of(4), 0u);
  BOOST_CHECK_EQUAL(table.size(), 2u);

  /* the bonds of a local particle */
  BOOST_CHECK_EQUAL(table.remove_bonds_of(2), 1u);
  BOOST_CHECK(table.view(1).empty());
  BOOST_CHECK_EQUAL(table.size(), 1u);

  /* a particle that left this rank is still known to its local partners */
  table.erase(3);
  BOOST_CHECK(table.empty());
  add(table, 0, {6, 7}, {6, 7});
  table.erase(7);
  BOOST_CHECK_EQUAL(table.size(), 1u);
  BOOST_CHECK_EQUAL(table.remove_bonds_of(7), 1u);
  BOOST_CHECK(table.empty());
  BOOST_CHECK_EQUAL(table.n_particles(), 0u);

  /* a particle that comes back is a local participant again */
  add(table, 0, {6, 7}, {6});
  std::vector<int> buffer{1, 0, 2, 6, 7};
  table.insert(7, buffer);
  BOOST_CHECK_EQUAL(table.size(), 1u);
  BOOST_CHECK_EQUAL(table.n_particles(), 2u);
  BOOST_CHECK((as_tuples(table.view(7)) == std::vector<Tuple>{{0, {6, 7}}}));
  BOOST_CHECK_EQUAL(table.remove_bonds_of(6), 1u);
  BOOST_CHECK(table.view(7).empty());
}

BOOST_AUTO_TEST_CASE(extract_and_insert) {
  /* rank A holds particles 1, 2, 3; rank B holds particle 4 */
  BondTable rank_a;
  BondTable rank_b;
  add(rank_a, 0, {1, 2}, {1, 2});
  add(rank_a, 1, {2, 1, 3}, {2, 1, 3});
  add(rank_a, 0, {3, 4}, {3});
  add(rank_b, 0, {3, 4}, {4});
  add(rank_b, 0, {4, 9}, {4});
  auto const bonds_of_2 = as_tuples(rank_a.view(2));
  auto const bonds_of_3 = as_tuples(rank_a.view(3));

  /* particles 3 and 2 move from A to B, in one message */
  std::vector<int> buffer;
  rank_a.extract(3, buffer);
  rank_a.extract(2, buffer);
  BOOST_CHECK(rank_a.view(2).empty());
  BOOST_CHECK(rank_a.view(3).empty());
  /* A keeps the bonds of particle 1 only */
  BOOST_CHECK((as_tuples(rank_a.view(1)) ==
               std::vector<Tuple>{{0, {1, 2}}, {1, {2, 1, 3}}}));
  BOOST_CHECK_EQUAL(rank_a.size(), 2u);

  auto rest = std::span<const int>(buffer);
  rest = rank_b.insert(3, rest);
  rest = rank_b.insert(2, rest);
  BOOST_CHECK(rest.empty());
  BOOST_CHECK((as_tuples(rank_b.view(2)) == bonds_of_2));
  BOOST_CHECK((as_tuples(rank_b.view(3)) == bonds_of_3));
  /* bonds already stored on B, or brought by both particles, are shared */
  BOOST_CHECK_EQUAL(rank_b.size(), 4u);
  BOOST_CHECK_EQUAL(n_live_rows(rank_b), 4u);
  BOOST_CHECK((as_tuples(rank_b.view(4)) ==
               std::vector<Tuple>{{0, {3, 4}}, {0, {4, 9}}}));

  /* a particle without bonds */
  std::vector<int> empty;
  rank_b.extract(42, empty);
  BondTable other;
  BOOST_CHECK(other.insert(42, empty).empty());
  BOOST_CHECK(other.empty());

  /* a copy does not remove the bonds */
  std::vector<int> copy;
  rank_b.copy_to(2, copy);
  BOOST_CHECK_EQUAL(rank_b.view(2).size(), 2u);
  other.insert(2, copy);
  BOOST_CHECK((as_tuples(other.view(2)) == bonds_of_2));
}

BOOST_AUTO_TEST_CASE(many_bonds) {
  /* the center of a star polymer, with all arms on this rank */
  BondTable table;
  int const n_arms = 40000;
  for (int i = 1; i <= n_arms; ++i) {
    table.add(0, std::vector<int>{0, i}, std::vector<int>{0, i});
  }
  BOOST_CHECK_EQUAL(table.size(), static_cast<std::size_t>(n_arms));
  BOOST_CHECK_EQUAL(table.view(0).size(), static_cast<std::size_t>(n_arms));
  BOOST_CHECK_EQUAL(table.view(17).size(), 1u);
  BOOST_CHECK(table.contains(0, std::vector<int>{0, n_arms}));
  BOOST_CHECK(table.remove(0, std::vector<int>{0, 1}));
  BOOST_CHECK_EQUAL(table.size(), static_cast<std::size_t>(n_arms - 1));
  BOOST_CHECK(table.view(1).empty());
}

BOOST_AUTO_TEST_CASE(malformed_input) {
  BondTable table;
  auto const check_throws = [&](std::vector<int> const &buffer) {
    BOOST_CHECK_THROW(table.insert(1, buffer), std::invalid_argument);
  };
  check_throws({});
  check_throws({-1});
  check_throws({1, 0, 2, 1});             // too short
  check_throws({1, 0, 5, 1, 2, 3, 4, 5}); // too many participants
  check_throws({1, 0, 0, 1, 2});          // no participants
  check_throws({1, -2, 2, 1, 2});         // invalid bond id
  check_throws({2, 0, 2, 1, 2, 0, 2, 1}); // second bond truncated
  /* a valid buffer is accepted */
  table.insert(1, std::vector<int>{1, 0, 2, 1, 2});
  BOOST_CHECK_EQUAL(table.size(), 1u);
}

BOOST_AUTO_TEST_CASE(inconsistent_state) {
  /* particle 2 no longer lists the bond (1, 2), but has other bonds */
  BondTable table;
  add(table, 0, {1, 2}, {1, 2});
  add(table, 0, {2, 3}, {2, 3});
  table.erase(2);
  add(table, 1, {2, 3}, {2, 3});
  BOOST_CHECK(remove(table, 0, {1, 2}));
  BOOST_CHECK(table.view(1).empty());
  BOOST_CHECK((as_tuples(table.view(2)) == std::vector<Tuple>{{1, {2, 3}}}));
  BOOST_CHECK((as_tuples(table.view(3)) ==
               std::vector<Tuple>{{0, {2, 3}}, {1, {2, 3}}}));
}

BOOST_AUTO_TEST_CASE(compaction) {
  BondTable table;
  int const n = 5000;
  for (int i = 0; i < n; ++i) {
    add(table, 0, {i, i + 1}, {i, i + 1});
  }
  /* remove most bonds */
  for (int i = 0; i < n; ++i) {
    if (i % 10 != 0) {
      BOOST_CHECK(remove(table, 0, {i, i + 1}));
    }
  }
  BOOST_CHECK(table.is_wasteful());
  auto const before = as_tuples(table.view(10));
  table.compact();
  BOOST_CHECK(not table.is_wasteful());
  BOOST_CHECK_EQUAL(table.rows().size(), table.size());
  BOOST_CHECK_EQUAL(table.size(), static_cast<std::size_t>(n / 10));
  /* the bonds are still found from every participant */
  BOOST_CHECK((as_tuples(table.view(10)) == before));
  BOOST_CHECK(table.contains(0, std::vector<int>{4990, 4991}));
  BOOST_CHECK(remove(table, 0, {4990, 4991}));
  BOOST_CHECK(table.view(4991).empty());
  add(table, 0, {7, 8}, {7, 8});
  BOOST_CHECK_EQUAL(table.size(), static_cast<std::size_t>(n / 10));
}
