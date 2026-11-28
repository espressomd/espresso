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

#include "BondTable.hpp"
#include "Particle.hpp"
#include "bond_queries.hpp"
#include "bonds.hpp"
#include "cell_system/CellStructure.hpp"
#include "cell_system/CellStructureType.hpp"
#include "communication.hpp"
#include "particle_node.hpp"
#include "system/System.hpp"

#include <utils/Vector.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <functional>
#include <iterator>
#include <ranges>
#include <vector>

/**
 * @brief Unit tests for add_bond()/remove_bond() (bonds.cpp) and the bond
 * queries (bond_queries.hpp).
 *
 * The bookkeeping of the bonds of one particle is covered by the
 * @ref BondTable unit tests; the interaction with MPI by the
 * collision-detection, bond-breakage, particle-slice and distributed-bond
 * regression tests.
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
auto const any_bond = [](BondRecord const &, int) { return true; };

/** @brief All bonds the particle @p pid participates in. */
auto bonds_of(int pid) {
  auto const pids = std::array<int, 1>{{pid}};
  return ::gather_bonds_of_particles(System::get_system(), pids, any_bond);
}

/** @brief Sorted ids of the bonds a particle owns. */
auto owned_bond_ids(int pid) {
  auto const pids = std::array<int, 1>{{pid}};
  std::vector<int> bond_ids;
  for (auto const &[bond_id, ids] : ::gather_bonds_of_particles(
           System::get_system(), pids, [](BondRecord const &bond, int p) {
             return bond.owner_id() == p;
           })) {
    bond_ids.push_back(bond_id);
  }
  std::ranges::sort(bond_ids);
  return bond_ids;
}

/** @brief Sorted ids of the bonds a particle participates in. */
auto all_bond_ids(int pid) {
  std::vector<int> bond_ids;
  for (auto const &[bond_id, ids] : bonds_of(pid)) {
    bond_ids.push_back(bond_id);
  }
  std::ranges::sort(bond_ids);
  return bond_ids;
}

/** @brief Total number of bonds, counted on their owners. */
std::size_t n_bonds() {
  std::size_t n_local = 0u;
  for (auto const &p : System::get_system().cell_structure->local_particles()) {
    n_local += static_cast<std::size_t>(std::ranges::distance(
        System::get_system().cell_structure->bond_table().owned(p.id())));
  }
  return boost::mpi::all_reduce(::comm_cart, n_local, std::plus<>());
}

/** @brief The particle, if it resides on this rank. */
Particle const *local(int pid) {
  auto const *p = System::get_system().cell_structure->get_local_particle(pid);
  return (p != nullptr and not p->is_ghost()) ? p : nullptr;
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
  BOOST_CHECK_EQUAL(n_bonds(), 0u);

  ::add_bond(system, pair_bond_id, {1, 2});
  ::add_bond(system, angle_bond_id, {2, 1, 3});
  BOOST_CHECK_EQUAL(n_bonds(), 2u);
  BOOST_CHECK(::is_bond_storage_consistent(system));

  /* every participant lists the bond, with the owner first */
  BOOST_CHECK(
      (bonds_of(1) == std::vector<BondTuple>{{pair_bond_id, {1, 2}},
                                             {angle_bond_id, {2, 1, 3}}}));
  BOOST_CHECK(
      (bonds_of(3) == std::vector<BondTuple>{{angle_bond_id, {2, 1, 3}}}));
  if (auto const *p = local(3)) {
    BOOST_CHECK_EQUAL(system.cell_structure->bond_table().view(p->id()).size(),
                      1u);
  }

  /* a bond can be removed in exact order, for all participants */
  ::remove_bond(system, angle_bond_id, {2, 1, 3});
  BOOST_CHECK_EQUAL(n_bonds(), 1u);
  BOOST_CHECK(bonds_of(3).empty());
  BOOST_CHECK((bonds_of(2) == std::vector<BondTuple>{{pair_bond_id, {1, 2}}}));
  ::remove_bond(system, pair_bond_id, {1, 2});
  BOOST_CHECK_EQUAL(n_bonds(), 0u);
  BOOST_CHECK(bonds_of(1).empty());
  BOOST_CHECK(bonds_of(2).empty());

  /* removing a bond that does not exist does nothing */
  ::add_bond(system, angle_bond_id, {2, 1, 3});
  ::remove_bond(system, pair_bond_id, {1, 2});
  ::remove_bond(system, angle_bond_id, {2, 3, 1});
  ::remove_bond(system, pair_bond_id, std::vector<int>{});
  BOOST_CHECK_EQUAL(n_bonds(), 1u);
  ::remove_bond(system, angle_bond_id, {2, 1, 3});
  BOOST_CHECK_EQUAL(n_bonds(), 0u);

  /* the same bond added from both sides: each removal removes in
   * exact order */
  ::add_bond(system, pair_bond_id, {1, 2});
  ::add_bond(system, pair_bond_id, {2, 1});
  BOOST_CHECK(::is_bond_storage_consistent(system));
  BOOST_CHECK_EQUAL(n_bonds(), 2u);
  ::remove_bond(system, pair_bond_id, {2, 1});
  BOOST_CHECK((bonds_of(2) == std::vector<BondTuple>{{pair_bond_id, {1, 2}}}));
  ::remove_bond(system, pair_bond_id, {1, 2});
  BOOST_CHECK_EQUAL(n_bonds(), 0u);

  /* a participant that no longer lists the bond is detected */
  ::add_bond(system, pair_bond_id, {1, 2});
  if (local(2)) {
    system.cell_structure->bond_table().erase(2);
  }
  BOOST_CHECK(not ::is_bond_storage_consistent(system));
}

BOOST_FIXTURE_TEST_CASE(bond_queries, RemoveAllParticles) {
  auto const pair_bond_id = 0;
  auto const angle_bond_id = 1;
  for (int pid = 1; pid <= 4; ++pid) {
    ::make_new_particle(pid, {1., 1., 1.});
  }
  auto &system = System::get_system();
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

  auto const *p1 = local(1);
  auto const *p2 = local(2);
  auto const *p3 = local(3);
  auto const *p4 = local(4);
  auto const &store = system.cell_structure->bond_table();
  if (p1 and p2 and p3 and p4) {
    /* the partners of an owned bond are in creation order */
    for (auto const bond : store.owned(2)) {
      BOOST_REQUIRE_EQUAL(bond.partner_ids().size(), 2u);
      BOOST_CHECK_EQUAL(bond.partner_ids()[0], 1);
      BOOST_CHECK_EQUAL(bond.partner_ids()[1], 3);
    }

    /* find_owned_bond_id() and find_bond_id_with() */
    auto const is_angle = [=](int id) { return id == angle_bond_id; };
    BOOST_CHECK(not find_owned_bond_id(store, *p1, is_angle));
    BOOST_CHECK(find_bond_id_with(store, *p1, is_angle) == angle_bond_id);
    BOOST_CHECK(find_owned_bond_id(store, *p2, is_angle) == angle_bond_id);
    BOOST_CHECK(not find_bond_id_with(store, *p4, is_angle));

    /* pair_bond_exists() finds a pair bond whichever particle owns it */
    BOOST_CHECK(pair_bond_exists(store, *p1, *p2, pair_bond_id));
    BOOST_CHECK(pair_bond_exists(store, *p2, *p1, pair_bond_id));
    BOOST_CHECK(not pair_bond_exists(store, *p1, *p3, pair_bond_id));
    BOOST_CHECK(not pair_bond_exists(store, *p1, *p1, pair_bond_id));
    /* angle bonds are not pair bonds */
    BOOST_CHECK(not pair_bond_exists(store, *p1, *p2, angle_bond_id));
    BOOST_CHECK(not pair_bond_exists(store, *p2, *p3, angle_bond_id));
    /* predicate overload */
    BOOST_CHECK(
        pair_bond_exists(store, *p2, *p1, [](int id) { return id >= 0; }));
    BOOST_CHECK(
        not pair_bond_exists(store, *p1, *p2, [](int) { return false; }));

    /* ghosts store no bonds: the bonds are read from the other particle */
    auto ghost = *p1;
    ghost.set_ghost(true);
    BOOST_CHECK(pair_bond_exists(store, ghost, *p2, pair_bond_id));
    BOOST_CHECK(pair_bond_exists(store, *p2, ghost, pair_bond_id));
  }
}

BOOST_FIXTURE_TEST_CASE(particle_removal, RemoveAllParticles) {
  auto const pair_bond_id = 0;
  auto const angle_bond_id = 1;
  for (int pid = 1; pid <= 4; ++pid) {
    ::make_new_particle(pid, {1., 1., 1.});
  }
  auto &system = System::get_system();
  ::add_bond(system, pair_bond_id, {1, 2});
  ::add_bond(system, angle_bond_id, {2, 1, 3});
  ::add_bond(system, pair_bond_id, {3, 4});
  BOOST_CHECK_EQUAL(n_bonds(), 3u);

  /* bonds of particles, restricted to the ones they own */
  auto const pids = std::array<int, 2>{{3, 4}};
  BOOST_CHECK((::gather_bonds_of_particles(
                   system, pids, [](BondRecord const &bond, int pid) {
                     return bond.owner_id() == pid;
                   }) == std::vector<BondTuple>{{pair_bond_id, {3, 4}}}));

  /* removing a particle removes all the bonds it participates in,
   * whether it owns them or not, on every rank */
  ::remove_particle(1);
  BOOST_CHECK(bonds_of(2).empty());
  BOOST_CHECK((bonds_of(3) == std::vector<BondTuple>{{pair_bond_id, {3, 4}}}));
  BOOST_CHECK_EQUAL(n_bonds(), 1u);
  BOOST_CHECK(::is_bond_storage_consistent(system));

  /* removing all particles removes all bonds */
  system.cell_structure->remove_all_particles();
  BOOST_CHECK_EQUAL(n_bonds(), 0u);
}

BOOST_AUTO_TEST_SUITE_END()
