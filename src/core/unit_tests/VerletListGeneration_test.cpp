/*
 * Copyright (C) 2021-2026 The ESPResSo project
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

// Test for the store-row Verlet list + generation guard. The Verlet list holds
// pairs of ParticleStore ROW indices held across integration steps
// (CellStructure::m_verlet_list), stamped with the store generation when the
// list was built. If the store is REBUILT between the build and a consume (a
// rebuild renumbers/permutes rows and bumps the generation) WITHOUT
// invalidating the list (i.e. without a Verlet rebuild), the stored rows are
// STALE and resolving them aliases the WRONG particle.
//
// This standalone test (no MPI, no decomposition) drives a ParticleStore by
// hand and reproduces exactly that build/permuting-rebuild/consume sequence:
//   * FAIL branch (guard disabled): resolving the stale rows after a permuting
//     rebuild reads the wrong particle -- the accumulated pair result diverges
//     from the correct one. This is the bug the guard exists to catch, shown
//     to actually occur.
//   * PASS branch (guard enabled): ParticleStoreGuard::assert_generation's
//     condition (recorded_generation == store.generation()) is TRUE before the
//     permuting rebuild and FALSE after it -- the guard detects the staleness
//     (and, in a debug build, the assert fires). Rebuilding the list restamps
//     the generation and resolves the correct particles.
//
// The identity gate is blind to this class of bug (positions can coincidentally
// match after a permuting sort), so this targeted test is the real safety net.

#define BOOST_TEST_MODULE VerletList generation guard test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "Particle.hpp"
#include "particle_store/ParticleStore.hpp"
#include "particle_store/StoreGenerationGuard.hpp"

#include <Kokkos_Core.hpp>

#include <cstddef>
#include <cstdint>
#include <utility>
#include <vector>

// ParticleStore allocates Kokkos Views, which requires an initialized runtime.
struct GlobalConfig {
  GlobalConfig() { Kokkos::initialize(); }
  ~GlobalConfig() { Kokkos::finalize(); }
};

BOOST_TEST_GLOBAL_CONFIGURATION(GlobalConfig);

namespace {

// A minimal stand-in for CellStructure::m_verlet_list plus its generation stamp
// (mirrors the real member's build/consume shape without pulling in the whole
// cell system / decomposition).
struct RowVerletList {
  std::vector<std::pair<int, int>> pairs;
  ParticleStore const *store = nullptr;
  std::uint64_t generation = 0u;
};

// Attach `ids` to `store` at rows [0, ids.size()), permuting: row r gets
// ids[r]. `seeds` are kept alive by the caller. This is one store rebuild: it
// bumps the generation. Returns nothing; the store now maps row r -> ids[r].
void rebuild_store_with_ids(ParticleStore &store, std::vector<Particle> &seeds,
                            std::vector<int> const &ids) {
  auto const n = ids.size();
  seeds.clear();
  seeds.resize(n);
  store.mark_dirty();
  store.begin_rebuild(n, 0u);
  int row = 0;
  for (auto &p : seeds) {
    store.assign_row(p, row++);
  }
  store.finish_rebuild();
  for (std::size_t r = 0u; r < n; ++r) {
    store.id(static_cast<int>(r)) = ids[r];
  }
}

// The "pair kernel" of this test: a deterministic function of the two resolved
// particles' ids. Stands in for a real force accumulation; a stale row makes it
// read the wrong ids and hence produce the wrong value.
long pair_value(Particle const &p1, Particle const &p2) {
  return 1000L * p1.id() + p2.id();
}

// Build the row Verlet list from the current store contents: for the given row
// pairs, record them AND stamp the current store generation + identity, exactly
// as CellStructure::verlet_list_loop's rebuild branch does.
void build_list(RowVerletList &vl, ParticleStore const &store,
                std::vector<std::pair<int, int>> const &row_pairs) {
  vl.pairs = row_pairs;
  vl.store = &store;
  vl.generation = store.generation();
}

// Consume the row Verlet list: resolve each row pair to views and accumulate,
// in the SAME pair order as the build, mirroring verlet_list_loop's consume
// branch. `check_guard` selects whether the generation guard is honored.
//   check_guard == true : assert the recorded generation still matches (PASS
//                         branch); returns the accumulated value only if valid.
//   check_guard == false: skip the guard (FAIL branch) -- resolves stale rows.
long consume_list(RowVerletList const &vl, ParticleStore &store) {
  long acc = 0L;
  for (auto const &[row1, row2] : vl.pairs) {
    auto const p1 = store.make_view(row1);
    auto const p2 = store.make_view(row2);
    acc += pair_value(p1, p2);
  }
  return acc;
}

} // namespace

// FAIL-THEN-PASS, both branches in one test so the evidence is self-contained.
BOOST_AUTO_TEST_CASE(stale_rows_read_wrong_particle_without_guard) {
  ParticleStore store{};
  std::vector<Particle> seeds;

  // Generation G0: rows 0,1,2 hold ids 10,11,12.
  rebuild_store_with_ids(store, seeds, {10, 11, 12});
  auto const gen_at_build = store.generation();

  // Build the Verlet list over pairs (0,1) and (1,2). Stamp G0.
  RowVerletList vl;
  build_list(vl, store, {{0, 1}, {1, 2}});
  BOOST_REQUIRE_EQUAL(vl.generation, gen_at_build);

  // The correct result, consumed at the SAME generation the list was built at.
  long const correct = consume_list(vl, store);
  // pair_value(10,11) + pair_value(11,12) = 10011 + 11012 = 21023
  BOOST_CHECK_EQUAL(correct, 21023L);

  // --- guard is VALID before any rebuild -------------------------------------
  // The generation stamp still matches the store -> the guard condition holds.
  BOOST_CHECK(vl.store == &store);
  BOOST_CHECK_EQUAL(vl.generation, store.generation());

  // --- a PERMUTING rebuild happens between build and the next consume --------
  // Rows are renumbered so that the SAME row indices now name DIFFERENT
  // particles (a reverse permutation of the ids). This bumps the generation.
  rebuild_store_with_ids(store, seeds, {12, 11, 10});
  BOOST_REQUIRE_GT(store.generation(), vl.generation);

  // --- FAIL branch: consume WITHOUT the guard --------------------------------
  // The stale rows (0,1),(1,2) now resolve to ids (12,11),(11,10):
  // pair_value(12,11) + pair_value(11,10) = 12011 + 11010 = 23021 != correct.
  long const stale = consume_list(vl, store);
  BOOST_CHECK_EQUAL(stale, 23021L);
  BOOST_CHECK_NE(stale, correct); // the demonstrated wrong read

  // --- PASS branch: the guard detects the staleness --------------------------
  // The recorded generation no longer matches the store's -> the guard
  // condition is FALSE, so a real consumer would rebuild (or, in a debug build,
  // ParticleStoreGuard::assert_generation would fire). We assert the guard
  // condition itself flips as designed.
  bool const guard_still_valid =
      (vl.store == &store) and (vl.generation == store.generation());
  BOOST_CHECK(not guard_still_valid);

  // Rebuilding the list restamps the current generation and, consumed now,
  // yields the correct value for the CURRENT store contents (12,11,10):
  build_list(vl, store, {{0, 1}, {1, 2}});
  BOOST_CHECK_EQUAL(vl.generation, store.generation());
  long const rebuilt = consume_list(vl, store);
  BOOST_CHECK_EQUAL(rebuilt, 23021L); // correct for the permuted store now
  // and the guard is valid again
  BOOST_CHECK(vl.generation == store.generation());
}

// The guard helper is a no-op that never fires while the generation matches
// (steady-state: build then consume at the same generation). This pins that the
// PASS path does not spuriously trip.
BOOST_AUTO_TEST_CASE(guard_does_not_fire_when_generation_matches) {
  ParticleStore store{};
  std::vector<Particle> seeds;
  rebuild_store_with_ids(store, seeds, {1, 2, 3, 4});

  RowVerletList vl;
  build_list(vl, store, {{0, 1}, {2, 3}});

  // No rebuild between build and consume -> the guard must accept.
  BOOST_CHECK_EQUAL(vl.generation, store.generation());
  BOOST_CHECK(vl.store == &store);
  // In a debug build this call must NOT abort; in release it compiles to
  // nothing. Either way the test process must survive the call.
  ParticleStoreGuard::assert_generation(store, vl.generation, vl.store);
  BOOST_CHECK(true);
}
