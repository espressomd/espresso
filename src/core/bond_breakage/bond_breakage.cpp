/*
 * Copyright (C) 2022 The ESPResSo project
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

#include "bond_breakage/bond_breakage.hpp"
#include "bond_breakage/actions.hpp"

#include "cell_system/CellStructure.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "system/System.hpp"

#include <utils/mpi/gather_buffer.hpp>
#include <utils/serialization/optional.hpp>

#include <boost/mpi.hpp>
#include <boost/serialization/access.hpp>
#include <boost/variant.hpp>

#include <cassert>
#include <memory>
#include <unordered_set>
#include <utility>
#include <vector>

namespace BondBreakage {

// Variant holding any of the actions
using Action = boost::variant<DeleteBond, DeleteAngleBond, DeleteAllBonds>;

// Set of actions
using ActionSet = std::unordered_set<Action>;

/** Add a particle+bond combination to the breakage queue */
void BondBreakage::queue_breakage(int particle_id,
                                  BondPartners const &bond_partners,
                                  int bond_type) {
  m_queue.emplace_back(QueueEntry{particle_id, bond_partners, bond_type});
}

/** @brief Gathers combined queue from all mpi ranks */
static auto gather_global_queue(Queue const &local_queue) {
  Queue res = local_queue;
  if (comm_cart.size() > 1) {
    Utils::Mpi::gather_buffer(res, comm_cart);
    boost::mpi::broadcast(comm_cart, res, 0);
  }
  return res;
}

/** @brief Constructs the actions to take for a breakage queue entry */
static ActionSet actions_for_breakage(CellStructure const &cell_structure,
                                      QueueEntry const &e,
                                      BreakageSpec const &spec) {
  auto is_angle_bond = [](auto const &bond_partners) {
    return bond_partners[1];
  }; // optional for second partner engaged

  // Handle different action types
  if (spec.action_type == ActionType::DELETE_BOND) {
    if (is_angle_bond(e.bond_partners)) {
      return {DeleteAngleBond{e.particle_id,
                              {{*(e.bond_partners[0]), *(e.bond_partners[1])}},
                              e.bond_type}};
    }
    return {DeleteBond{e.particle_id, *(e.bond_partners[0]), e.bond_type}};
  }
#ifdef VIRTUAL_SITES_RELATIVE
  // revert bind at point of collision for pair bonds
  if (spec.action_type == ActionType::REVERT_BIND_AT_POINT_OF_COLLISION and
      not is_angle_bond(e.bond_partners)) {
    // We need to find the base particles for the two virtual sites
    // between which the bond broke.
    auto p1 = cell_structure.get_local_particle(e.particle_id);
    auto p2 = cell_structure.get_local_particle(*(e.bond_partners[0]));
    if (p1 and p2) {
      if (not p1->is_virtual() or not p2->is_virtual()) {
        runtimeErrorMsg() << "The REVERT_BIND_AT_POINT_OF_COLLISION bond "
                             "breakage action has to be configured for the "
                             "bond on the virtual site. Encountered a particle "
                             "that is not virtual.";
        return {};
      }

      return {
          // Bond between virtual sites
          DeleteBond{e.particle_id, *(e.bond_partners[0]), e.bond_type},
          // Bond between base particles. We do not know, on which of these
          // the bond is defined, since bonds are stored only on one partner
          DeleteAllBonds{p1->vs_relative().to_particle_id,
                         p2->vs_relative().to_particle_id},
          DeleteAllBonds{p2->vs_relative().to_particle_id,
                         p1->vs_relative().to_particle_id},
      };
    }
  } else {
    // revert bind at point of collision for angle bonds
    auto vs = cell_structure.get_local_particle(e.particle_id);
    auto p1 = cell_structure.get_local_particle(*(e.bond_partners[0]));
    auto p2 = cell_structure.get_local_particle(*(e.bond_partners[1]));
    if (p1 and p2) {
      if (not vs->is_virtual()) {
        runtimeErrorMsg() << "The REVERT_BIND_AT_POINT_OF_COLLISION bond "
                             "breakage action has to be configured for the "
                             "bond on the virtual site. Encountered a particle "
                             "that is not virtual.";
        return {};
      }

      return {// Angle bond on the virtual site
              DeleteAngleBond{e.particle_id, {p1->id(), p2->id()}, e.bond_type},
              // Bond between base particles. We do not know, on which of these
              // the bond is defined, since bonds are stored only on one partner
              DeleteAllBonds{p1->id(), p2->id()},
              DeleteAllBonds{p2->id(), p1->id()}};
    }
  }
#endif // VIRTUAL_SITES_RELATIVE
  return {};
}

/**
 * @brief Delete specific bond.
 */
static void remove_bond(Particle &p, BondView const &view) {
  auto &bond_list = p.bonds();
  auto it = std::find(bond_list.begin(), bond_list.end(), view);
  if (it != bond_list.end()) {
    bond_list.erase(it);
  }
}

/**
 * @brief Delete pair bonds to a specific partner
 */
static void remove_pair_bonds_to(Particle &p, int other_pid) {
  std::vector<std::pair<int, int>> to_delete;
  for (auto b : p.bonds()) {
    if (b.partner_ids().size() == 1 and b.partner_ids()[0] == other_pid)
      to_delete.emplace_back(b.bond_id(), other_pid);
  }
  for (auto const &b : to_delete) {
    remove_bond(p, BondView(b.first, {&b.second, 1}));
  }
}

// Handler for the different delete events
class execute : public boost::static_visitor<> {
  CellStructure &cell_structure;

public:
  explicit execute(CellStructure &cell_structure)
      : cell_structure{cell_structure} {}

  void operator()(DeleteBond const &d) const {
    if (auto p = cell_structure.get_local_particle(d.particle_id)) {
      remove_bond(*p, BondView(d.bond_type, {&d.bond_partner_id, 1}));
    }
  }
  void operator()(DeleteAngleBond const &d) const {
    if (auto p = cell_structure.get_local_particle(d.particle_id)) {
      remove_bond(*p, BondView(d.bond_type, {&d.bond_partner_id[0], 2}));
    }
  }
  void operator()(DeleteAllBonds const &d) const {
    if (auto p = cell_structure.get_local_particle(d.particle_id_1)) {
      remove_pair_bonds_to(*p, d.particle_id_2);
    }
  }
};

void BondBreakage::process_queue_impl(System::System &system) {
  auto global_queue = gather_global_queue(m_queue);
  auto &cell_structure = *system.cell_structure;

  // Construct delete actions from breakage queue
  ActionSet actions = {};
  for (auto const &e : global_queue) {
    // Retrieve relevant breakage spec
    assert(breakage_specs.count(e.bond_type) != 0);
    auto const &spec = breakage_specs.at(e.bond_type);
    actions.merge(actions_for_breakage(cell_structure, e, *spec));
  }

  // Execute actions
  for (auto const &a : actions) {
    boost::apply_visitor(execute(cell_structure), a);
    system.on_particle_change();
  }
}

bool bond_handler(BondBreakage &bond_breakage, Particle &p,
                  std::span<Particle *> partners, int bond_id,
                  BoxGeometry const &box_geo) {

  if (partners.size() == 1u) { // pair bonds
    auto d = box_geo.get_mi_vector(p.pos(), partners[0]->pos()).norm();
    if (bond_breakage.check_and_handle_breakage(
            p.id(), {{partners[0]->id(), std::nullopt}}, bond_id, d)) {
      return true;
    }
    return false;
  }
  if (partners.size() == 2u) { // angle bond
    auto d =
        box_geo.get_mi_vector(partners[0]->pos(), partners[1]->pos()).norm();
    if (bond_breakage.check_and_handle_breakage(
            p.id(), {{partners[0]->id(), partners[1]->id()}}, bond_id, d)) {
      return true;
    }
    return false;
  }
  return false;
}

void execute_bond_breakage(System::System &system,
                           BondBreakage &bond_breakage) {
  // Clear the bond breakage queue
  bond_breakage.clear_queue();

  // Create the bond kernel function (the bond handler)
  auto bond_kernel = [&](Particle &p, int bond_id,
                         std::span<Particle *> partners) {
    return bond_handler(bond_breakage, p, partners, bond_id, *system.box_geo);
  };

  // Use the CellStructure::bond_loop to process bonds
  system.cell_structure->bond_loop(bond_kernel);

  // Process the bond breakage queue
  bond_breakage.process_queue(system);
}

} // namespace BondBreakage
