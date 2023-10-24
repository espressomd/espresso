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

#pragma once

#include "system/System.hpp"

#include <boost/serialization/access.hpp>

#include <array>
#include <memory>
#include <optional>
#include <unordered_map>
#include <vector>

namespace BondBreakage {

/** Stores one or two bond parnters for pair/angle bonds */
using BondPartners = std::array<std::optional<int>, 2>;

enum class ActionType {
  NONE = 0,
  DELETE_BOND = 1,
  REVERT_BIND_AT_POINT_OF_COLLISION = 2
};

struct BreakageSpec {
  double breakage_length;
  ActionType action_type;
};

// Broken bond record
struct QueueEntry {
  int particle_id;
  BondPartners bond_partners = {};
  int bond_type;

  // Serialization for synchronization across mpi ranks
  friend class boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar &particle_id;
    ar &bond_partners;
    ar &bond_type;
  }
};

/** @brief Record bonds broken during a time step. */
using Queue = std::vector<QueueEntry>;

class BondBreakage {
  Queue m_queue;

public:
  /** @brief Bond breakage specifications. */
  std::unordered_map<int, std::shared_ptr<BreakageSpec>> breakage_specs;

  BondBreakage() : m_queue{}, breakage_specs{} {}

  /** @brief Check if the bond between the particles should break, if yes, queue
   *  it.
   */
  bool check_and_handle_breakage(int particle_id,
                                 BondPartners const &bond_partners,
                                 int bond_type, double distance) {
    if (breakage_specs.count(bond_type) == 0) {
      return false; // No breakage rule for this bond type
    }

    // Retrieve relevant breakage spec
    auto const &spec = *(breakage_specs.at(bond_type));

    // Is the bond length longer than the breakage length?
    if (distance >= spec.breakage_length) {
      queue_breakage(particle_id, bond_partners, bond_type);
      return true;
    }
    return false;
  }

  void clear_queue() { m_queue.clear(); }

  void process_queue(System::System &system) {
    if (not breakage_specs.empty()) {
      process_queue_impl(system);
    }
  }

private:
  void process_queue_impl(System::System &system);

  /** Add a particle+bond combination to the breakage queue */
  void queue_breakage(int particle_id, BondPartners const &bond_partners,
                      int bond_type);
};

} // namespace BondBreakage
