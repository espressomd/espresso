/*
 * Copyright (C) 2025-2026 The ESPResSo project
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

#define BOOST_TEST_MODULE "bonds"
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "EspressoCoreGlobalConfig.hpp"
#include "ParticleFactory.hpp"

#include "BondList.hpp"
#include "bonds.hpp"
#include "cell_system/CellStructure.hpp"
#include "cell_system/CellStructureType.hpp"
#include "system/System.hpp"

#include <utils/Vector.hpp>

#include <algorithm>
#include <array>
#include <vector>

/**
 * @brief Unit tests for add_bond()/remove_bond() (bonds.cpp).
 *
 * These exercise the "primary entry + mirror entries on every other
 * participant" bookkeeping and, in particular, remove_bond()'s two-phase
 * matching strategy (try particle_ids[0] as the primary owner first, else
 * fall back to matching any role independently per participant): a single
 * bond can be removed starting from any of its participants, not just the
 * one holding the primary entry, and this must not misbehave when the same
 * bond id/participant set exists twice with swapped ownership. All
 * particles here are co-located, so a single MPI rank sees every
 * participant as a real (non-ghost) local particle; the interaction of this
 * bookkeeping with ghost communication across ranks is covered separately
 * by the collision-detection and particle-slice-assignment regression
 * tests (which specifically pin bond partners on different ranks).
 */

struct GlobalConfig : public EspressoCoreGlobalConfig {
  GlobalConfig() {
    auto system = System::System::create();
    System::set_system(system);
    system->set_cell_structure_topology(CellStructureType::REGULAR);
    system->set_box_l(Utils::Vector3d::broadcast(10.));
  }
  ~GlobalConfig() { System::reset_system(); }
};

BOOST_TEST_GLOBAL_CONFIGURATION(GlobalConfig);
BOOST_AUTO_TEST_SUITE(suite)

namespace {
auto sorted_partner_ids(BondView const &bond) {
  auto ids =
      std::vector<int>(bond.partner_ids().begin(), bond.partner_ids().end());
  std::ranges::sort(ids);
  return ids;
}

/** @brief All bond entries currently on the local particle @p pid,
 *  or an empty vector if it does not exist locally. */
auto local_bonds(int pid) {
  auto &system = System::get_system();
  auto bonds = std::vector<BondView>{};
  if (auto const *p = system.cell_structure->get_local_particle(pid)) {
    for (auto const bond : p->bonds()) {
      bonds.push_back(bond);
    }
  }
  return bonds;
}
} // namespace

BOOST_FIXTURE_TEST_CASE(remove_from_owner_side, ParticleFactory) {
  auto const bond_id = 0;
  create_particle({1., 1., 1.}, 1, 0);
  create_particle({1., 1., 1.}, 2, 0);

  auto &system = System::get_system();
  ::add_bond(system, bond_id, {1, 2});

  BOOST_REQUIRE_EQUAL(local_bonds(1).size(), 1u);
  BOOST_REQUIRE_EQUAL(local_bonds(2).size(), 1u);
  BOOST_CHECK(local_bonds(1).front().is_primary());
  BOOST_CHECK(not local_bonds(2).front().is_primary());

  /* remove_bond() called with the owner (1) first removes both the
   * primary entry on 1 and the mirror entry on 2. */
  BOOST_CHECK(::remove_bond(system, bond_id, {1, 2}));
  BOOST_CHECK(local_bonds(1).empty());
  BOOST_CHECK(local_bonds(2).empty());
}

BOOST_FIXTURE_TEST_CASE(remove_from_mirror_side, ParticleFactory) {
  auto const bond_id = 0;
  create_particle({1., 1., 1.}, 1, 0);
  create_particle({1., 1., 1.}, 2, 0);

  auto &system = System::get_system();
  ::add_bond(system, bond_id, {1, 2});
  BOOST_REQUIRE_EQUAL(local_bonds(1).size(), 1u);
  BOOST_REQUIRE_EQUAL(local_bonds(2).size(), 1u);

  /* remove_bond() called with the mirror holder (2) first must still
   * find and remove both entries: particle_ids[0]=2 does not own a
   * primary entry for this bond, so the primary-first fast path fails
   * and the fallback ("match any role independently per participant")
   * must kick in for both 2's own mirror and 1's primary. */
  BOOST_CHECK(::remove_bond(system, bond_id, {2, 1}));
  BOOST_CHECK(local_bonds(1).empty());
  BOOST_CHECK(local_bonds(2).empty());
}

BOOST_FIXTURE_TEST_CASE(remove_via_mirror_only_participant_of_angle_bond,
                        ParticleFactory) {
  /* A 3-participant (angle) bond owned by 1, with 2 and 3 only ever
   * holding mirrors -- neither 2 nor 3 has ever been party to a primary
   * entry for any bond, unlike the pair-bond case above. */
  auto const bond_id = 0;
  create_particle({1., 1., 1.}, 1, 0);
  create_particle({1., 1., 1.}, 2, 0);
  create_particle({1., 1., 1.}, 3, 0);

  auto &system = System::get_system();
  ::add_bond(system, bond_id, {1, 2, 3});

  BOOST_REQUIRE_EQUAL(local_bonds(1).size(), 1u);
  BOOST_REQUIRE_EQUAL(local_bonds(2).size(), 1u);
  BOOST_REQUIRE_EQUAL(local_bonds(3).size(), 1u);
  BOOST_CHECK(local_bonds(1).front().is_primary());
  BOOST_CHECK(not local_bonds(2).front().is_primary());
  BOOST_CHECK(not local_bonds(3).front().is_primary());
  BOOST_CHECK(
      (sorted_partner_ids(local_bonds(2).front()) == std::vector<int>{1, 3}));
  BOOST_CHECK(
      (sorted_partner_ids(local_bonds(3).front()) == std::vector<int>{1, 2}));

  /* Initiate removal from 3, a participant that only ever held a
   * mirror; the primary owner (1) is listed afterwards. */
  BOOST_CHECK(::remove_bond(system, bond_id, {3, 1, 2}));
  BOOST_CHECK(local_bonds(1).empty());
  BOOST_CHECK(local_bonds(2).empty());
  BOOST_CHECK(local_bonds(3).empty());
}

BOOST_FIXTURE_TEST_CASE(duplicate_bond_added_from_both_sides, ParticleFactory) {
  /* The same bond id/participant pair added once from each side creates
   * two co-existing entries per particle (a primary and a mirror). A
   * single remove_bond() call naming an owner must remove exactly the
   * matching primary/mirror pair, not an arbitrary (potentially
   * mismatched) one -- verified by checking what is left behind, not
   * just that something was removed. */
  auto const bond_id = 0;
  create_particle({1., 1., 1.}, 1, 0);
  create_particle({1., 1., 1.}, 2, 0);

  auto &system = System::get_system();
  ::add_bond(system, bond_id, {1, 2}); // "bond A": primary on 1, mirror on 2
  ::add_bond(system, bond_id, {2, 1}); // "bond B": primary on 2, mirror on 1

  BOOST_REQUIRE_EQUAL(local_bonds(1).size(), 2u);
  BOOST_REQUIRE_EQUAL(local_bonds(2).size(), 2u);

  /* Remove bond A by naming its owner (1) first: must remove 1's
   * primary and 2's matching mirror, leaving bond B (2's primary, 1's
   * mirror) completely untouched on both particles. */
  BOOST_CHECK(::remove_bond(system, bond_id, {1, 2}));

  BOOST_REQUIRE_EQUAL(local_bonds(1).size(), 1u);
  BOOST_REQUIRE_EQUAL(local_bonds(2).size(), 1u);
  BOOST_CHECK(not local_bonds(1).front().is_primary());
  BOOST_CHECK(local_bonds(2).front().is_primary());

  /* Remove the remaining bond B by naming its owner (2) first. */
  BOOST_CHECK(::remove_bond(system, bond_id, {2, 1}));
  BOOST_CHECK(local_bonds(1).empty());
  BOOST_CHECK(local_bonds(2).empty());
}

BOOST_FIXTURE_TEST_CASE(rebuild_bond_mirrors_backfills_missing_mirrors,
                        ParticleFactory) {
  /* rebuild_bond_mirrors() exists to backfill mirrors for checkpoints
   * written before bonds were stored on all participants, which only
   * ever contain primary entries (mpiio.cpp calls it, unconditionally,
   * right after loading particle data). Simulate that on-disk state
   * directly: insert primary-only entries via BondList::insert(),
   * bypassing add_bond() so no mirrors are written, exactly like
   * loading such a legacy archive into p.bonds() would leave things. */
  auto const pair_bond_id = 0;
  auto const angle_bond_id = 1;
  create_particle({1., 1., 1.}, 1, 0);
  create_particle({1., 1., 1.}, 2, 0);
  create_particle({1., 1., 1.}, 3, 0);

  auto &system = System::get_system();
  auto &cell_structure = *system.cell_structure;

  {
    auto const partner = std::array<int, 1>{{2}};
    cell_structure.get_local_particle(1)->bonds().insert(
        BondView{pair_bond_id, partner, true});
  }
  {
    auto const partners = std::array<int, 2>{{1, 3}};
    cell_structure.get_local_particle(2)->bonds().insert(
        BondView{angle_bond_id, partners, true});
  }

  BOOST_REQUIRE_EQUAL(local_bonds(1).size(), 1u); // pair primary
  BOOST_REQUIRE_EQUAL(local_bonds(2).size(), 1u); // angle primary
  BOOST_CHECK(local_bonds(3).empty());            // no mirror yet

  ::rebuild_bond_mirrors(system);

  /* The pair bond's mirror now exists on 2, and the angle bond's
   * mirrors now exist on 1 and 3, matching what add_bond() would have
   * produced for the same participant lists. */
  BOOST_REQUIRE_EQUAL(local_bonds(1).size(), 2u); // pair primary + angle mirror
  BOOST_REQUIRE_EQUAL(local_bonds(2).size(), 2u); // angle primary + pair mirror
  BOOST_REQUIRE_EQUAL(local_bonds(3).size(), 1u); // angle mirror

  auto const has_bond = [](std::vector<BondView> const &bonds, int bond_id,
                           bool is_primary, std::vector<int> const &partners) {
    return std::ranges::any_of(bonds, [&](BondView const &b) {
      return b.bond_id() == bond_id and b.is_primary() == is_primary and
             sorted_partner_ids(b) == partners;
    });
  };
  BOOST_CHECK(has_bond(local_bonds(1), pair_bond_id, true, {2}));
  BOOST_CHECK(has_bond(local_bonds(1), angle_bond_id, false, {2, 3}));
  BOOST_CHECK(has_bond(local_bonds(2), angle_bond_id, true, {1, 3}));
  BOOST_CHECK(has_bond(local_bonds(2), pair_bond_id, false, {1}));
  BOOST_CHECK(has_bond(local_bonds(3), angle_bond_id, false, {1, 2}));

  /* Existing mirrors are left untouched: calling it again must not
   * duplicate anything. */
  ::rebuild_bond_mirrors(system);
  BOOST_CHECK_EQUAL(local_bonds(1).size(), 2u);
  BOOST_CHECK_EQUAL(local_bonds(2).size(), 2u);
  BOOST_CHECK_EQUAL(local_bonds(3).size(), 1u);
}

BOOST_AUTO_TEST_SUITE_END()
