/*
 * Copyright (C) 2017-2026 The ESPResSo project
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

// Tests for the range-collapsed cell iteration surface: a Cell holds a
// contiguous (offset, count) store-row range (CellRowSpan), and
// RowParticleRange yields Particle views over it -- skipping rows marked
// pending-removed on the store. All are standalone (no MPI, no decomposition):
// a ParticleStore is built by hand and a Cell's range is set the way
// ensure_particle_store_synchronized's write-back does, so the tests pin the
// machinery's contract without the integration path (which the physics ctest
// battery exercises end-to-end under ADDITIONAL_CHECKS).

#define BOOST_TEST_MODULE RowParticleRange test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "cell_system/Cell.hpp"
#include "cell_system/CellRows.hpp"
#include "cell_system/RowParticleRange.hpp"

#include "Particle.hpp"
#include "particle_store/ParticleStore.hpp"

#include <Kokkos_Core.hpp>

#include <cstddef>
#include <vector>

// ParticleStore allocates Kokkos Views, which requires an initialized runtime.
struct GlobalConfig {
  GlobalConfig() { Kokkos::initialize(); }
  ~GlobalConfig() { Kokkos::finalize(); }
};

BOOST_TEST_GLOBAL_CONFIGURATION(GlobalConfig);

namespace {
// Build a store with `n` rows and set row r's id to `first_id + r`. The seed
// particles must be kept alive by the caller (they back the store columns);
// the caller only needs the store afterwards.
void build_store(ParticleStore &store, std::vector<Particle> &seeds,
                 std::size_t n_local, std::size_t n_ghost, int first_id) {
  auto const n = n_local + n_ghost;
  seeds.clear();
  seeds.resize(n);
  store.mark_dirty();
  store.begin_rebuild(n_local, n_ghost);
  int row = 0;
  for (auto &p : seeds) {
    store.assign_row(p, row++);
  }
  store.finish_rebuild();
  for (int r = 0; r < static_cast<int>(n); ++r) {
    store.id(r) = first_id + r;
  }
}

// Wire a cell to the store and give it the contiguous store-row range
// [offset, offset+count), the way ensure_particle_store_synchronized's
// (offset, count) write-back does.
void fill_cell(Cell &cell, ParticleStore &store, std::size_t offset,
               std::size_t count) {
  cell.set_store(store);
  cell.set_range(offset, count);
}
} // namespace

// The cell's committed range records the exact store rows, in ascending order.
// Mirrors the ensure_particle_store_synchronized write-back for a "local" and a
// "ghost" cell block.
BOOST_AUTO_TEST_CASE(cell_range_matches_store_order) {
  ParticleStore store{};
  std::vector<Particle> seeds;
  // 3 local rows [0,3) then 2 ghost rows [3,5), ids 100..104.
  build_store(store, seeds, 3u, 2u, 100);

  Cell local_cell, ghost_cell;
  fill_cell(local_cell, store, 0u, 3u);
  fill_cell(ghost_cell, store, 3u, 2u);

  auto check = [](Cell const &cell, int base_id) {
    auto const bag = cell.particles();
    auto const rows = cell.rows();
    BOOST_REQUIRE_EQUAL(rows.size(), bag.size());
    int expected = base_id;
    for (auto const &p : bag) {
      // store.id(row) == view.id(), the ADDITIONAL_CHECKS invariant.
      BOOST_CHECK_EQUAL(p.id(), expected);
      ++expected;
    }
  };
  check(local_cell, 100);
  check(ghost_cell, 103);

  BOOST_REQUIRE_EQUAL(local_cell.rows().size(), 3u);
  BOOST_CHECK_EQUAL(local_cell.offset(), 0u);
  BOOST_CHECK_EQUAL(local_cell.count(), 3u);
  BOOST_REQUIRE_EQUAL(ghost_cell.rows().size(), 2u);
  BOOST_CHECK_EQUAL(ghost_cell.offset(), 3u);
  BOOST_CHECK_EQUAL(ghost_cell.count(), 2u);
}

// The range yields the cell's particles by ascending row, ids matching the
// store, on a built (multi-particle) store.
BOOST_AUTO_TEST_CASE(row_range_multi) {
  ParticleStore store{};
  std::vector<Particle> seeds;
  build_store(store, seeds, 4u, 0u, 200);

  Cell cell;
  fill_cell(cell, store, 0u, 4u);

  std::vector<int> range_ids;
  RowParticleRange range{cell.rows(), store};
  BOOST_CHECK_EQUAL(range.size(), cell.rows().size());
  for (auto const &p : range) {
    range_ids.push_back(p.id());
  }

  BOOST_REQUIRE_EQUAL(range_ids.size(), 4u);
  BOOST_CHECK_EQUAL(range_ids[0], 200);
  BOOST_CHECK_EQUAL(range_ids[3], 203);
}

// A SINGLE-row cell yields exactly one particle.
BOOST_AUTO_TEST_CASE(row_range_single) {
  ParticleStore store{};
  std::vector<Particle> seeds;
  build_store(store, seeds, 2u, 0u, 400);

  Cell cell;
  fill_cell(cell, store, 1u, 1u); // just row 1 (id 401)
  RowParticleRange range{cell.rows(), store};
  BOOST_CHECK_EQUAL(range.size(), 1u);
  std::size_t count = 0u;
  for (auto const &p : range) {
    BOOST_CHECK_EQUAL(p.id(), 401);
    ++count;
  }
  BOOST_CHECK_EQUAL(count, 1u);
}

// An EMPTY range yields an empty range (begin() == end()).
BOOST_AUTO_TEST_CASE(row_range_empty) {
  ParticleStore store{};
  std::vector<Particle> seeds;
  build_store(store, seeds, 1u, 0u, 0);

  Cell cell; // range defaults to (0, 0)
  cell.set_store(store);
  RowParticleRange range{cell.rows(), store};
  BOOST_CHECK_EQUAL(range.size(), 0u);
  BOOST_CHECK(range.begin() == range.end());
  std::size_t count = 0u;
  for (auto const &p : range) {
    static_cast<void>(p);
    ++count;
  }
  BOOST_CHECK_EQUAL(count, 0u);
}

// PENDING-REMOVAL visibility semantics: a row marked pending-removed is
// invisible to iteration and to the live size, BEFORE any rebuild resolves it.
// Marking the middle / first / last row each drops exactly that particle from
// the range while the range keeps its store-row order for the survivors.
BOOST_AUTO_TEST_CASE(row_range_skips_pending_removed) {
  ParticleStore store{};
  std::vector<Particle> seeds;
  build_store(store, seeds, 5u, 0u, 500); // ids 500..504

  Cell cell;
  fill_cell(cell, store, 0u, 5u);
  BOOST_REQUIRE_EQUAL(cell.rows().size(), 5u);
  BOOST_REQUIRE_EQUAL(cell.size(), 5u);

  // Drop the middle row (row 2, id 502).
  store.mark_pending_removal(2);
  BOOST_CHECK_EQUAL(cell.rows().size(), 4u); // live count excludes it
  BOOST_CHECK_EQUAL(cell.size(), 4u);        // Cell::size() excludes it too
  {
    std::vector<int> ids;
    for (auto const &p : cell.particles()) {
      ids.push_back(p.id());
    }
    BOOST_REQUIRE_EQUAL(ids.size(), 4u);
    BOOST_CHECK_EQUAL(ids[0], 500);
    BOOST_CHECK_EQUAL(ids[1], 501);
    BOOST_CHECK_EQUAL(ids[2], 503); // 502 skipped, order preserved
    BOOST_CHECK_EQUAL(ids[3], 504);
  }
  // The CellRowSpan iteration yields the raw surviving rows, skipping row 2.
  {
    std::vector<int> rows;
    for (int const r : cell.rows()) {
      rows.push_back(r);
    }
    std::vector<int> const expected{0, 1, 3, 4};
    BOOST_CHECK(rows == expected);
  }

  // Also drop the first (row 0) and last (row 4): only 501, 503 remain.
  store.mark_pending_removal(0);
  store.mark_pending_removal(4);
  {
    std::vector<int> ids;
    for (auto const &p : cell.particles()) {
      ids.push_back(p.id());
    }
    std::vector<int> const expected{501, 503};
    BOOST_CHECK(ids == expected);
  }

  // Dropping every row yields an empty range.
  store.mark_pending_removal(1);
  store.mark_pending_removal(3);
  {
    RowParticleRange range{cell.rows(), store};
    BOOST_CHECK_EQUAL(range.size(), 0u);
    BOOST_CHECK(range.begin() == range.end());
  }
}

// LIFETIME CONTRACT: `Particle &p = *it` stays valid while the iterator lives
// and is not incremented. Binding a reference to the view, then reading it
// after other statements, must still refer to the same particle. This is the
// pattern the bond handlers rely on (Particle &p2 = *ptr held across a call).
BOOST_AUTO_TEST_CASE(
    row_range_dereference_reference_stable_while_iterator_lives) {
  ParticleStore store{};
  std::vector<Particle> seeds;
  build_store(store, seeds, 3u, 0u, 300);

  Cell cell;
  fill_cell(cell, store, 0u, 3u);

  RowParticleRange range{cell.rows(), store};
  auto it = range.begin();

  // Bind a reference to the current view and keep using it across other
  // statements without re-dereferencing the iterator.
  Particle &p = *it;
  int const id_at_bind = p.id();
  BOOST_CHECK_EQUAL(id_at_bind, 300);

  // Do unrelated work; the reference must remain valid and refer to the same
  // particle (row 0), i.e. the referent was not a destroyed temporary.
  volatile int scratch = 0;
  for (int i = 0; i < 1000; ++i) {
    // Note: not `scratch += i;` -- compound assignment to a volatile operand
    // is deprecated in C++20 and rejected under -Werror by gcc-12.
    scratch = scratch + i;
  }
  static_cast<void>(scratch);

  BOOST_CHECK_EQUAL(p.id(), 300);
  BOOST_CHECK_EQUAL(p.store_row(), 0);
  BOOST_CHECK(p.store() == &store);
  // A write through the still-live reference reaches the store row it names.
  p.pos() = Utils::Vector3d{1., 2., 3.};
  BOOST_CHECK_EQUAL(store.position_value(0)[1], 2.);

  // After incrementing, the SAME iterator now names the next row; a fresh
  // dereference reflects the advance (the cache is refreshed, per the
  // documented contract that incrementing invalidates the previous referent).
  ++it;
  Particle &p_next = *it;
  BOOST_CHECK_EQUAL(p_next.id(), 301);
  BOOST_CHECK_EQUAL(p_next.store_row(), 1);
}
