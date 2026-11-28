/*
 * Copyright (C) 2022-2026 The ESPResSo project
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

#define BOOST_TEST_MODULE "Bond breakage"
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "config/config.hpp"

#include "EspressoCoreGlobalConfig.hpp"

#include "bond_breakage/actions.hpp"
#include "bond_breakage/bond_breakage.hpp"

#include "Particle.hpp"
#include "cell_system/CellStructure.hpp"
#include "cell_system/CellStructureType.hpp"
#include "system/System.hpp"

#include <memory>
#include <optional>
#include <utility>

struct GlobalConfig : public EspressoCoreGlobalConfig {
  ~GlobalConfig() { ::System::reset_system(); }
};

BOOST_TEST_GLOBAL_CONFIGURATION(GlobalConfig);

BOOST_AUTO_TEST_CASE(test_actions_equality) {
  {
    using Action = BondBreakage::DeleteBond;
    BOOST_CHECK((Action{1, 2, 3} == Action{1, 2, 3}));
    BOOST_CHECK((Action{1, 2, 3} != Action{0, 2, 3}));
    BOOST_CHECK((Action{1, 2, 3} != Action{1, 0, 3}));
    BOOST_CHECK((Action{1, 2, 3} != Action{1, 2, 0}));
  }

  {
    using Action = BondBreakage::DeleteAngleBond;
    BOOST_CHECK((Action{1, {2, 4}, 3} == Action{1, {2, 4}, 3}));
    BOOST_CHECK((Action{1, {2, 4}, 3} != Action{0, {2, 4}, 3}));
    BOOST_CHECK((Action{1, {2, 4}, 3} != Action{1, {0, 4}, 3}));
    BOOST_CHECK((Action{1, {2, 4}, 3} != Action{1, {2, 0}, 3}));
    BOOST_CHECK((Action{1, {2, 4}, 3} != Action{1, {2, 4}, 0}));
  }

  {
    using Action = BondBreakage::DeleteAllBonds;
    BOOST_CHECK((Action{1, 2} == Action{1, 2}));
    BOOST_CHECK((Action{1, 2} != Action{0, 2}));
    BOOST_CHECK((Action{1, 2} != Action{1, 0}));
  }
}

BOOST_AUTO_TEST_CASE(test_actions_hash_value) {
  {
    using Action = BondBreakage::DeleteBond;
    BOOST_CHECK((Action{1, 2, 3}.hash_value() == Action{1, 2, 3}.hash_value()));
    BOOST_CHECK((Action{1, 2, 3}.hash_value() != Action{0, 2, 3}.hash_value()));
    BOOST_CHECK((Action{1, 2, 3}.hash_value() != Action{1, 0, 3}.hash_value()));
    BOOST_CHECK((Action{1, 2, 3}.hash_value() != Action{1, 2, 0}.hash_value()));
  }

  {
    // clang-format off
    using Action = BondBreakage::DeleteAngleBond;
    BOOST_CHECK((Action{1, {2, 4}, 3}.hash_value() == Action{1, {2, 4}, 3}.hash_value()));
    BOOST_CHECK((Action{1, {2, 4}, 3}.hash_value() != Action{0, {2, 4}, 3}.hash_value()));
    BOOST_CHECK((Action{1, {2, 4}, 3}.hash_value() != Action{1, {0, 4}, 3}.hash_value()));
    BOOST_CHECK((Action{1, {2, 4}, 3}.hash_value() != Action{1, {2, 0}, 3}.hash_value()));
    BOOST_CHECK((Action{1, {2, 4}, 3}.hash_value() != Action{1, {2, 4}, 0}.hash_value()));
    // clang-format on
  }

  {
    using Action = BondBreakage::DeleteAllBonds;
    BOOST_CHECK((Action{1, 2}.hash_value() == Action{1, 2}.hash_value()));
    BOOST_CHECK((Action{1, 2}.hash_value() != Action{0, 2}.hash_value()));
    BOOST_CHECK((Action{1, 2}.hash_value() != Action{1, 0}.hash_value()));
  }
}

#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE
// Add a particle with the given id to the cell structure.
static void add_particle(CellStructure &cell_structure, int id) {
  Particle p;
  p.id() = id;
  cell_structure.add_particle(std::move(p));
}

/**
 * The breakage queue is broadcast to all ranks. On a rank which holds both
 * bond partners but not the central virtual site, the angle-bond branch of
 * @ref BondBreakage::actions_for_breakage must produce no action.
 * This test reproduces that locality condition on a single rank.
 */
BOOST_AUTO_TEST_CASE(actions_for_breakage_angle_missing_central_particle) {
  auto system = System::System::create();
  System::set_system(system);
  system->set_cell_structure_topology(CellStructureType::NSQUARE);
  auto &cell_structure = *system->cell_structure;

  auto const central_id = 0;
  auto const partner_id_1 = 1;
  auto const partner_id_2 = 2;

  // both partners are present, the central particle is absent
  add_particle(cell_structure, partner_id_1);
  add_particle(cell_structure, partner_id_2);

  BOOST_REQUIRE(cell_structure.get_local_particle(central_id) == nullptr);
  BOOST_REQUIRE(cell_structure.get_local_particle(partner_id_1) != nullptr);
  BOOST_REQUIRE(cell_structure.get_local_particle(partner_id_2) != nullptr);

  BondBreakage::QueueEntry entry{
      central_id,
      {{std::optional<int>{partner_id_1}, std::optional<int>{partner_id_2}}},
      /* bond_type */ 0};
  BondBreakage::BreakageSpec spec{
      /* breakage_length */ 0.,
      BondBreakage::ActionType::REVERT_BIND_AT_POINT_OF_COLLISION};

  auto const actions =
      BondBreakage::actions_for_breakage(cell_structure, entry, spec);
  BOOST_CHECK_EQUAL(actions.size(), 0u);

  System::reset_system();
}
#endif // ESPRESSO_VIRTUAL_SITES_RELATIVE
