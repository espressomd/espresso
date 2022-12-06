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

#include "cells.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "particle_data.hpp"

#include <utils/mpi/gather_buffer.hpp>

#include <boost/mpi.hpp>
#include <boost/optional.hpp>
#include <boost/serialization/access.hpp>
#include <boost/variant.hpp>

#include <cassert>
#include <memory>
#include <unordered_set>
#include <utility>
#include <vector>

namespace BondBreakage {

// Bond breakage specifications
static std::unordered_map<int, std::shared_ptr<BreakageSpec>> breakage_specs;

void insert_spec(int key, std::shared_ptr<BreakageSpec> obj) {
  breakage_specs[key] = std::move(obj);
}

void erase_spec(int key) { breakage_specs.erase(key); }

// Variant holding any of the actions
using Action = boost::variant<DeleteBond, DeleteAllBonds>;

// Set of actions
using ActionSet = std::unordered_set<Action>;

// Broken bond record
struct QueueEntry {
  int particle_id;
  int bond_partner_id;
  int bond_type;

  // Serialization for synchronization across mpi ranks
  friend class boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar &particle_id;
    ar &bond_partner_id;
    ar &bond_type;
  }
};

/** @brief Queue to record bonds broken during a time step */
using Queue = std::vector<QueueEntry>;
static Queue queue;

/** @brief Retrieve breakage specification for the bond type */
boost::optional<BreakageSpec> get_breakage_spec(int bond_type) {
  if (breakage_specs.find(bond_type) != breakage_specs.end()) {
    return {*(breakage_specs.at(bond_type))};
  }
  return {};
}

/** Add a particle+bond combination to the breakage queue */
void queue_breakage(int particle_id, int bond_partner_id, int bond_type) {
  queue.emplace_back(QueueEntry{particle_id, bond_partner_id, bond_type});
}

bool check_and_handle_breakage(int particle_id, int bond_partner_id,
                               int bond_type, double distance) {
  // Retrieve specification for this bond type
  auto spec = get_breakage_spec(bond_type);
  if (!spec)
    return false; // No breakage rule for this bond type

  // Is the bond length longer than the breakage length?
  if (distance >= (*spec).breakage_length) {
    queue_breakage(particle_id, bond_partner_id, bond_type);
    return true;
  }
  return false;
}

void clear_queue() { queue.clear(); }

/** @brief Gathers combined queue from all mpi ranks */
Queue gather_global_queue(Queue const &local_queue) {
  Queue res = local_queue;
  Utils::Mpi::gather_buffer(res, comm_cart);
  boost::mpi::broadcast(comm_cart, res, 0);
  return res;
}

/** @brief Constructs the actions to take for a breakage queue entry */
ActionSet actions_for_breakage(QueueEntry const &e) {
  // Retrieve relevant breakage spec
  auto const spec = get_breakage_spec(e.bond_type);
  assert(spec);

  // Handle different action types
  if ((*spec).action_type == ActionType::DELETE_BOND)
    return {DeleteBond{e.particle_id, e.bond_partner_id, e.bond_type}};
#ifdef VIRTUAL_SITES_RELATIVE
  if ((*spec).action_type == ActionType::REVERT_BIND_AT_POINT_OF_COLLISION) {
    // We need to find the base particles for the two virtual sites
    // between which the bond broke.
    auto p1 = cell_structure.get_local_particle(e.particle_id);
    auto p2 = cell_structure.get_local_particle(e.bond_partner_id);
    if (not p1 or not p2)
      return {}; // particles not on this mpi rank

    if (not p1->is_virtual() or not p2->is_virtual()) {
      runtimeErrorMsg() << "The REVERT_BIND_AT_POINT_OF_COLLISION bond "
                           "breakage action has to be configured for the "
                           "bond on the virtual site. Encountered a particle "
                           "that is not virtual.";
      return {};
    }

    return {
        // Bond between virtual sites
        DeleteBond{e.particle_id, e.bond_partner_id, e.bond_type},
        // Bond between base particles. We do not know, on which of the two
        // the bond is defined, since bonds are stored only on one partner
        DeleteAllBonds{p1->vs_relative().to_particle_id,
                       p2->vs_relative().to_particle_id},
        DeleteAllBonds{p2->vs_relative().to_particle_id,
                       p1->vs_relative().to_particle_id},
    };
  }
#endif // VIRTUAL_SITES_RELATIVE
  return {};
}

// Handler for the different delete events
class execute : public boost::static_visitor<> {
public:
  void operator()(DeleteBond const &d) const {
    auto p = cell_structure.get_local_particle(d.particle_id);
    if (!p)
      return;
    local_remove_bond(*p, {d.bond_type, d.bond_partner_id});
  }
  void operator()(DeleteAllBonds const &d) const {
    auto p = cell_structure.get_local_particle(d.particle_id_1);
    if (!p)
      return;
    local_remove_pair_bonds_to(*p, d.particle_id_2);
  }
};

void process_queue() {
  if (breakage_specs.empty())
    return;

  auto global_queue = gather_global_queue(queue);

  // Construct delete actions from breakage queue
  ActionSet actions = {};
  for (auto const &e : global_queue) {
    actions.merge(actions_for_breakage(e));
  }

  // Execute actions
  for (auto const &a : actions) {
    boost::apply_visitor(execute(), a);
  }
}
} // namespace BondBreakage
