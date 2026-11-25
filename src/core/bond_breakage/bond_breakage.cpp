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

#include <config/config.hpp>

#include "bond_breakage/actions.hpp"
#include "bond_breakage/bond_breakage.hpp"

#include "bonds.hpp"
#include "cell_system/CellStructure.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "system/System.hpp"

#include <utils/mpi/gather_buffer.hpp>
#include <utils/serialization/optional.hpp>

#include <boost/mpi.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/std_variant.hpp>
#include <boost/serialization/vector.hpp>

#include <cassert>
#include <memory>
#include <mutex>
#include <span>
#include <unordered_set>
#include <utility>
#include <variant>
#include <vector>

namespace BondBreakage {

// Variant holding any of the actions
using Action = std::variant<DeleteBond, DeleteAngleBond, DeleteAllBonds>;

// Set of actions
using ActionSet = std::unordered_set<Action>;

/** Add a particle+bond combination to the breakage queue */
void BondBreakage::queue_breakage(int particle_id,
                                  BondPartners const &bond_partners,
                                  int bond_type) {
  {
    std::lock_guard<std::mutex> lock(queue_mtx);
    m_queue.emplace_back(QueueEntry{particle_id, bond_partners, bond_type});
  }
}

/** @brief Combine the elements collected on all MPI ranks, in rank order,
 *  into the same vector on every rank. */
template <typename T>
static std::vector<T> gather_on_all_ranks(std::vector<T> elements) {
  if (comm_cart.size() > 1) {
    Utils::Mpi::gather_buffer(elements, comm_cart);
    boost::mpi::broadcast(comm_cart, elements, 0);
  }
  return elements;
}

/** @brief Constructs the actions to take for a breakage queue entry */
static ActionSet actions_for_breakage(CellStructure const &cell_structure,
                                      QueueEntry const &e,
                                      BreakageSpec const &spec) {
  auto is_angle_bond = [](auto const &bond_partners) {
    return bond_partners[1];
  }; // optional for second partner engaged

#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE
  if (spec.action_type == ActionType::REVERT_BIND_AT_POINT_OF_COLLISION) {
    if (not is_angle_bond(e.bond_partners)) {
      // revert bind at point of collision for pair bonds
      // We need to find the base particles for the two virtual sites
      // between which the bond broke.
      auto p1 = cell_structure.get_local_particle(e.particle_id);
      auto p2 = cell_structure.get_local_particle(*(e.bond_partners[0]));
      if (p1 and p2) {
        if (not p1->is_virtual() or not p2->is_virtual()) {
          runtimeErrorMsg()
              << "The REVERT_BIND_AT_POINT_OF_COLLISION bond "
                 "breakage action has to be configured for the "
                 "bond on the virtual site. Encountered a particle "
                 "that is not virtual.";
          return {};
        }

        return {
            // Bond between virtual sites
            DeleteBond{e.particle_id, *(e.bond_partners[0]), e.bond_type},
            // Bond between base particles, of unspecified bond type.
            DeleteAllBonds{p1->vs_relative().to_particle_id,
                           p2->vs_relative().to_particle_id},
        };
      }
    } else {
      // revert bind at point of collision for angle bonds
      auto vs = cell_structure.get_local_particle(e.particle_id);
      auto p1 = cell_structure.get_local_particle(*(e.bond_partners[0]));
      auto p2 = cell_structure.get_local_particle(*(e.bond_partners[1]));
      if (p1 and p2) {
        if (not vs->is_virtual()) {
          runtimeErrorMsg()
              << "The REVERT_BIND_AT_POINT_OF_COLLISION bond "
                 "breakage action has to be configured for the "
                 "bond on the virtual site. Encountered a particle "
                 "that is not virtual.";
          return {};
        }

        return {
            // Angle bond on the virtual site
            DeleteAngleBond{e.particle_id, {p1->id(), p2->id()}, e.bond_type},
            // Bond between base particles, of unspecified bond type.
            DeleteAllBonds{p1->id(), p2->id()}};
      }
    }
  }
#endif // ESPRESSO_VIRTUAL_SITES_RELATIVE

  assert(spec.action_type == ActionType::DELETE_BOND);
  if (is_angle_bond(e.bond_partners)) {
    return {DeleteAngleBond{e.particle_id,
                            {{*(e.bond_partners[0]), *(e.bond_partners[1])}},
                            e.bond_type}};
  }
  return {DeleteBond{e.particle_id, *(e.bond_partners[0]), e.bond_type}};
}

// Handler for the different delete events
class execute {
  System::System &system;

public:
  explicit execute(System::System &system) : system{system} {}

  void operator()(DeleteBond const &d) const {
    ::remove_bond(system, d.bond_type, {d.particle_id, d.bond_partner_id});
  }
  void operator()(DeleteAngleBond const &d) const {
    ::remove_bond(system, d.bond_type,
                  {d.particle_id, d.bond_partner_id[0], d.bond_partner_id[1]});
  }
  void operator()(DeleteAllBonds const &d) const {
    // Remove every pair bond (any type, every occurrence) between the two
    // particles, owned by either one, from the bond table of this rank, if
    // one of the two particles resides here. Executed identically on all
    // ranks, so that the bonds are removed on every rank that stores them.
    auto const id1 = d.particle_id_1;
    auto const id2 = d.particle_id_2;
    auto const connects = [id1, id2](BondRecord const &bond) {
      auto const ids = bond.particle_ids();
      return ids.size() == 2u and ((ids[0] == id1 and ids[1] == id2) or
                                   (ids[0] == id2 and ids[1] == id1));
    };
    for (auto const pid : {id1, id2}) {
      auto *p = system.cell_structure->get_local_particle(pid);
      if (p != nullptr and not p->is_ghost()) {
        system.cell_structure->bond_table().remove_if(pid, connects);
      }
    }
  }
};

void BondBreakage::process_queue_impl(System::System &system) {
  auto const global_queue = gather_on_all_ranks(m_queue);
  auto &cell_structure = *system.cell_structure;

  if (global_queue.empty()) {
    return;
  }

  // Construct delete actions from breakage queue
  ActionSet local_actions = {};
  for (auto const &e : global_queue) {
    // Retrieve relevant breakage spec
    assert(breakage_specs.contains(e.bond_type));
    auto const &spec = breakage_specs.at(e.bond_type);
    local_actions.merge(actions_for_breakage(cell_structure, e, *spec));
  }

  // Some actions can only be constructed on the ranks where the involved
  // particles are known, but the bonds must be removed on every rank where
  // a participant resides: execute the actions found on any rank, on every
  // rank and in the same order.
  auto const actions = gather_on_all_ranks(
      std::vector<Action>(local_actions.begin(), local_actions.end()));
  ActionSet executed = {};
  for (auto const &action : actions) {
    if (executed.insert(action).second) {
      std::visit(execute(system), action);
    }
  }
  if (not executed.empty()) {
    system.on_particle_change();
  }
}

static bool bond_handler(BondBreakage &bond_breakage, Particle &p,
                         std::span<Particle *> partners, int bond_id,
                         BoxGeometry const &box_geo) {
  auto retval = false;
  if (partners.size() == 1u) { // pair bonds
    auto d = box_geo.get_mi_vector(p.pos(), partners[0]->pos()).norm();
    retval = bond_breakage.check_and_handle_breakage(
        p.id(), {{partners[0]->id(), std::nullopt}}, bond_id, d);
  } else if (partners.size() == 2u) { // angle bond
    auto d =
        box_geo.get_mi_vector(partners[0]->pos(), partners[1]->pos()).norm();
    retval = bond_breakage.check_and_handle_breakage(
        p.id(), {{partners[0]->id(), partners[1]->id()}}, bond_id, d);
  }
  return retval;
}

void BondBreakage::execute_bond_breakage(System::System &system) {
  system.cell_structure->update_ghosts_and_resort_particle(
      system.get_global_ghost_flags());

  // Clear the bond breakage queue
  clear_queue();

  // Create the bond kernel function (the bond handler)
  auto bond_kernel = [&](Particle &p, int bond_id,
                         std::span<Particle *> partners) {
    bond_handler(*this, p, partners, bond_id, *system.box_geo);
    return false;
  };

  // Use the CellStructure::bond_loop to process bonds
  system.cell_structure->bond_loop(bond_kernel);

  // Process the bond breakage queue
  process_queue(system);
}

} // namespace BondBreakage
