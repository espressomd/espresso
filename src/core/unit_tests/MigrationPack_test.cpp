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

// The maximally-populated round-trip in this file is the FIELD-COMPLETENESS
// ENFORCEMENT for the per-field migration pack and for
// ParticleStore::copy_row. Every ifdef-guarded store field is written to a
// distinct sentinel value in fill_maximal() and compared field-for-field in
// check_row_equal(); a field missed by pack/unpack or by copy_row FAILS a test
// rather than being silently dropped.

#define BOOST_TEST_MODULE MigrationPack test
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <config/config.hpp>

#include "ParticleStoreMaximalPopulation.hpp"

#include "particle_store/MigrationPack.hpp"
#include "particle_store/ParticleStore.hpp"

#include "Particle.hpp"

#include <Kokkos_Core.hpp>

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

// ParticleStore allocates Kokkos Views, which requires an initialized runtime.
struct GlobalConfig {
  GlobalConfig() { Kokkos::initialize(); }
  ~GlobalConfig() { Kokkos::finalize(); }
};

BOOST_TEST_GLOBAL_CONFIGURATION(GlobalConfig);

// The maximal-population field-completeness helpers (make_store,
// fill_maximal, check_row_equal) live in a shared header so the permute
// rebuild tests reuse the same enforcement.
using namespace maximal_population;

// (a) pack -> unpack round-trip identity on a maximally-populated store.
// Destination store starts empty/different-size; verify field-for-field
// equality incl. ragged contents. THIS is the field-completeness enforcement.
BOOST_AUTO_TEST_CASE(pack_unpack_roundtrip_maximal) {
  constexpr std::size_t n = 3u;
  auto source = make_store(n);
  for (std::size_t r = 0u; r < n; ++r) {
    fill_maximal(source, static_cast<int>(r), 1000. + 7. * double(r));
  }

  std::array<int, n> const rows{0, 1, 2};
  std::vector<char> buffer;
  MigrationPack::pack_rows(source, rows, buffer);

  // Destination starts DIFFERENT size (5 rows) and all-default; only the first
  // n rows are overwritten by unpack.
  auto dest = make_store(5u);
  auto const count = MigrationPack::unpack_rows(dest, 0, buffer);
  BOOST_CHECK_EQUAL(count, n);

  for (std::size_t r = 0u; r < n; ++r) {
    check_row_equal(source, static_cast<int>(r), dest, static_cast<int>(r));
  }
}

// unpack into a non-zero first_row (models a receiver appending migrated rows
// after its existing locals).
BOOST_AUTO_TEST_CASE(pack_unpack_roundtrip_offset_first_row) {
  auto source = make_store(2u);
  fill_maximal(source, 0, 2000.);
  fill_maximal(source, 1, 2100.);

  std::array<int, 2> const rows{0, 1};
  std::vector<char> buffer;
  MigrationPack::pack_rows(source, rows, buffer);

  auto dest = make_store(5u);
  auto const count = MigrationPack::unpack_rows(dest, 3, buffer);
  BOOST_CHECK_EQUAL(count, 2u);
  check_row_equal(source, 0, dest, 3);
  check_row_equal(source, 1, dest, 4);
}

// (b) copy_row completeness — same maximal row copied between stores; verify
// every field. Independent path from pack/unpack, same enforcement helper.
BOOST_AUTO_TEST_CASE(copy_row_completeness_across_stores) {
  auto source = make_store(2u);
  fill_maximal(source, 0, 3000.);
  fill_maximal(source, 1, 3100.);

  auto dest = make_store(4u);
  dest.copy_row(source, 0, 2);
  dest.copy_row(source, 1, 3);

  check_row_equal(source, 0, dest, 2);
  check_row_equal(source, 1, dest, 3);
}

// copy_row within a single store (self-copy: source == *this).
BOOST_AUTO_TEST_CASE(copy_row_within_one_store) {
  auto store = make_store(3u);
  fill_maximal(store, 0, 4000.);
  store.copy_row(store, 0, 2);
  check_row_equal(store, 0, store, 2);
}

// (c) edge case: zero rows -> only the row-count header (u64, value 0).
BOOST_AUTO_TEST_CASE(pack_unpack_zero_rows) {
  auto source = make_store(1u);
  fill_maximal(source, 0, 5000.);

  std::array<int, 0> const rows{};
  std::vector<char> buffer;
  MigrationPack::pack_rows(source, rows, buffer);
  BOOST_CHECK_EQUAL(buffer.size(), sizeof(std::uint64_t)); // just the header

  auto dest = make_store(1u);
  auto const count = MigrationPack::unpack_rows(dest, 0, buffer);
  BOOST_CHECK_EQUAL(count, 0u);
}

// (c) edge case: bond-free particles mixed with bonded ones (and, where
// compiled, exclusion-free mixed with exclusion-bearing).
BOOST_AUTO_TEST_CASE(pack_unpack_mixed_ragged_presence) {
  constexpr std::size_t n = 4u;
  auto source = make_store(n);
  // row 0: bonds + exclusions; row 1: neither; row 2: bonds only;
  // row 3: exclusions only.
  fill_maximal(source, 0, 6000., /*with_bonds=*/true, /*with_exclusions=*/true);
  fill_maximal(source, 1, 6100., /*with_bonds=*/false,
               /*with_exclusions=*/false);
  fill_maximal(source, 2, 6200., /*with_bonds=*/true,
               /*with_exclusions=*/false);
  fill_maximal(source, 3, 6300., /*with_bonds=*/false,
               /*with_exclusions=*/true);

  std::array<int, n> const rows{0, 1, 2, 3};
  std::vector<char> buffer;
  MigrationPack::pack_rows(source, rows, buffer);

  auto dest = make_store(n);
  auto const count = MigrationPack::unpack_rows(dest, 0, buffer);
  BOOST_CHECK_EQUAL(count, n);
  for (std::size_t r = 0u; r < n; ++r) {
    check_row_equal(source, static_cast<int>(r), dest, static_cast<int>(r));
  }
  // Spell out the bond-free row's emptiness explicitly.
  BOOST_CHECK_EQUAL(dest.bonds_sidecar_reference(1).size(), 0u);
#ifdef ESPRESSO_EXCLUSIONS
  BOOST_CHECK_EQUAL(dest.exclusions_sidecar_reference(1).size(), 0u);
#endif
}

// (d) packed_size == the actual buffer size pack_rows produces, both with and
// without ragged data present.
BOOST_AUTO_TEST_CASE(packed_size_matches_buffer) {
  auto source = make_store(3u);
  fill_maximal(source, 0, 7000., /*with_bonds=*/true, /*with_exclusions=*/true);
  fill_maximal(source, 1, 7100., /*with_bonds=*/false,
               /*with_exclusions=*/false);
  fill_maximal(source, 2, 7200., /*with_bonds=*/true, /*with_exclusions=*/true);

  std::array<int, 3> const rows{0, 1, 2};
  auto const predicted = MigrationPack::packed_size(source, rows);
  std::vector<char> buffer;
  MigrationPack::pack_rows(source, rows, buffer);
  BOOST_CHECK_EQUAL(predicted, buffer.size());

  // Empty selection: header only.
  std::array<int, 0> const none{};
  BOOST_CHECK_EQUAL(MigrationPack::packed_size(source, none),
                    sizeof(std::uint64_t));
}
