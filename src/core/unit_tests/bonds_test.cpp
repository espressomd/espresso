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

#include "BondStorage.hpp"
#include "bond_queries.hpp"
#include "bonds.hpp"
#include "cell_system/CellStructure.hpp"
#include "cell_system/CellStructureType.hpp"
#include "particle_node.hpp"
#include "system/System.hpp"

#include <utils/Vector.hpp>

#include <algorithm>
#include <vector>

/**
 * @brief Unit tests for add_bond()/remove_bond() (bonds.cpp) and the bond
 * queries (bond_queries.hpp).
 *
 * The bookkeeping of the bond topology itself is covered by the
 * @ref BondStorage unit tests; the interaction with MPI by the
 * collision-detection, bond-breakage and particle-slice regression tests.
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
auto const any_bond = [](BondStorage::Bond const &) { return true; };

/** @brief All bonds the particle @p pid participates in. */
auto bonds_of(int pid) {
  return ::get_bonds_of_particle(System::get_system(), pid, any_bond);
}

/** @brief Sorted ids of the bonds a particle owns. */
auto owned_bond_ids(int pid) {
  std::vector<int> bond_ids;
  System::get_system().bond_storage->for_each_owned_bond(
      pid,
      [&](BondStorage::Bond const &bond) { bond_ids.push_back(bond.bond_id); });
  std::ranges::sort(bond_ids);
  return bond_ids;
}

/** @brief Sorted ids of the bonds a particle participates in. */
auto all_bond_ids(int pid) {
  std::vector<int> bond_ids;
  System::get_system().bond_storage->for_each_bond_with(
      pid,
      [&](BondStorage::Bond const &bond) { bond_ids.push_back(bond.bond_id); });
  std::ranges::sort(bond_ids);
  return bond_ids;
}

struct RemoveAllParticles {
  ~RemoveAllParticles() {
    System::get_system().cell_structure->remove_all_particles();
  }
};
} // namespace

BOOST_FIXTURE_TEST_CASE(add_and_remove, RemoveAllParticles) {
  auto const pair_bond_id = 0;
  auto const angle_bond_id = 1;
  for (int pid = 1; pid <= 3; ++pid) {
    ::make_new_particle(pid, {1., 1., 1.});
  }
  auto &system = System::get_system();
  auto const &bond_storage = *system.bond_storage;
  BOOST_CHECK(bond_storage.empty());

  ::add_bond(system, pair_bond_id, {1, 2});
  ::add_bond(system, angle_bond_id, {2, 1, 3});
  BOOST_CHECK_EQUAL(bond_storage.size(), 2u);
  BOOST_CHECK(::is_bond_storage_consistent(system));

  /* a bond is found from every participant, with the owner first */
  BOOST_CHECK(
      (bonds_of(1) == std::vector<BondTuple>{{pair_bond_id, {1, 2}},
                                             {angle_bond_id, {2, 1, 3}}}));
  BOOST_CHECK(
      (bonds_of(3) == std::vector<BondTuple>{{angle_bond_id, {2, 1, 3}}}));

  /* a bond can be removed from any participant, in any order */
  BOOST_CHECK(::remove_bond(system, angle_bond_id, {3, 2, 1}));
  BOOST_CHECK(bonds_of(3).empty());
  BOOST_CHECK(::remove_bond(system, pair_bond_id, {2, 1}));
  BOOST_CHECK(bond_storage.empty());

  /* removing a bond that does not exist is reported */
  BOOST_CHECK(not ::remove_bond(system, pair_bond_id, {1, 2}));

  /* the same bond added from both sides: each removal removes one of
   * them, preferring the one owned by the first particle */
  ::add_bond(system, pair_bond_id, {1, 2});
  ::add_bond(system, pair_bond_id, {2, 1});
  BOOST_CHECK(::remove_bond(system, pair_bond_id, {2, 1}));
  BOOST_CHECK((bonds_of(2) == std::vector<BondTuple>{{pair_bond_id, {1, 2}}}));
  BOOST_CHECK(::remove_bond(system, pair_bond_id, {2, 1}));
  BOOST_CHECK(bond_storage.empty());
}

BOOST_FIXTURE_TEST_CASE(bond_queries, RemoveAllParticles) {
  auto const pair_bond_id = 0;
  auto const angle_bond_id = 1;
  for (int pid = 1; pid <= 3; ++pid) {
    ::make_new_particle(pid, {1., 1., 1.});
  }
  auto &system = System::get_system();
  auto const &bond_storage = *system.bond_storage;
  ::add_bond(system, pair_bond_id, {1, 2});
  ::add_bond(system, angle_bond_id, {2, 1, 3});

  /* owned bonds: every bond is visited exactly once over all particles,
   * on its owner */
  BOOST_CHECK((owned_bond_ids(1) == std::vector<int>{pair_bond_id}));
  BOOST_CHECK((owned_bond_ids(2) == std::vector<int>{angle_bond_id}));
  BOOST_CHECK(owned_bond_ids(3).empty());
  /* all bonds: a bond is visited once per participant */
  BOOST_CHECK(
      (all_bond_ids(1) == std::vector<int>{pair_bond_id, angle_bond_id}));
  BOOST_CHECK(
      (all_bond_ids(2) == std::vector<int>{pair_bond_id, angle_bond_id}));
  BOOST_CHECK((all_bond_ids(3) == std::vector<int>{angle_bond_id}));
  BOOST_CHECK(all_bond_ids(4).empty());

  /* the partners of an owned bond are in creation order */
  bond_storage.for_each_owned_bond(2, [&](BondStorage::Bond const &bond) {
    BOOST_REQUIRE_EQUAL(bond.partner_ids().size(), 2u);
    BOOST_CHECK_EQUAL(bond.partner_ids()[0], 1);
    BOOST_CHECK_EQUAL(bond.partner_ids()[1], 3);
  });

  /* find_owned_bond_id() and find_bond_id_with() */
  auto const is_angle = [=](int id) { return id == angle_bond_id; };
  BOOST_CHECK(not find_owned_bond_id(bond_storage, 1, is_angle));
  BOOST_CHECK(find_bond_id_with(bond_storage, 1, is_angle) == angle_bond_id);
  BOOST_CHECK(find_owned_bond_id(bond_storage, 2, is_angle) == angle_bond_id);
  BOOST_CHECK(not find_bond_id_with(bond_storage, 4, is_angle));

  /* pair_bond_exists() finds a pair bond whichever particle owns it */
  BOOST_CHECK(pair_bond_exists(bond_storage, 1, 2, pair_bond_id));
  BOOST_CHECK(pair_bond_exists(bond_storage, 2, 1, pair_bond_id));
  BOOST_CHECK(not pair_bond_exists(bond_storage, 1, 3, pair_bond_id));
  BOOST_CHECK(not pair_bond_exists(bond_storage, 1, 1, pair_bond_id));
  /* angle bonds are not pair bonds */
  BOOST_CHECK(not pair_bond_exists(bond_storage, 1, 2, angle_bond_id));
  BOOST_CHECK(not pair_bond_exists(bond_storage, 2, 3, angle_bond_id));
  /* predicate overload */
  BOOST_CHECK(
      pair_bond_exists(bond_storage, 2, 1, [](int id) { return id >= 0; }));
  BOOST_CHECK(
      not pair_bond_exists(bond_storage, 1, 2, [](int) { return false; }));
}

BOOST_FIXTURE_TEST_CASE(particle_removal, RemoveAllParticles) {
  auto const pair_bond_id = 0;
  auto const angle_bond_id = 1;
  for (int pid = 1; pid <= 4; ++pid) {
    ::make_new_particle(pid, {1., 1., 1.});
  }
  auto &system = System::get_system();
  auto const &bond_storage = *system.bond_storage;
  ::add_bond(system, pair_bond_id, {1, 2});
  ::add_bond(system, angle_bond_id, {2, 1, 3});
  ::add_bond(system, pair_bond_id, {3, 4});
  BOOST_CHECK_EQUAL(bond_storage.size(), 3u);

  /* bonds of a particle, restricted to the ones it owns */
  BOOST_CHECK(
      (::get_bonds_of_particle(system, 3, [](BondStorage::Bond const &bond) {
         return bond.owner_id() == 3;
       }) == std::vector<BondTuple>{{pair_bond_id, {3, 4}}}));

  /* removing a particle removes all the bonds it participates in,
   * whether it owns them or not */
  ::remove_particle(1);
  BOOST_CHECK(bonds_of(1).empty());
  BOOST_CHECK(bonds_of(2).empty());
  BOOST_CHECK((bonds_of(3) == std::vector<BondTuple>{{pair_bond_id, {3, 4}}}));
  BOOST_CHECK_EQUAL(bond_storage.size(), 1u);
  BOOST_CHECK(::is_bond_storage_consistent(system));

  /* removing all particles removes all bonds */
  system.cell_structure->remove_all_particles();
  BOOST_CHECK(bond_storage.empty());
}

BOOST_AUTO_TEST_SUITE_END()
