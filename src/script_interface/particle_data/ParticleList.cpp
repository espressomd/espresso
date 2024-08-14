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

#include "ParticleList.hpp"
#include "ParticleHandle.hpp"
#include "ParticleSlice.hpp"

#include "script_interface/ObjectState.hpp"
#include "script_interface/ScriptInterface.hpp"
#include "script_interface/system/Leaf.hpp"

#include "core/cell_system/CellStructure.hpp"
#include "core/exclusions.hpp"
#include "core/particle_node.hpp"
#include "core/system/System.hpp"

#include <utils/Vector.hpp>
#include <utils/mpi/gather_buffer.hpp>
#include <utils/serialization/pack.hpp>

#include <boost/mpi/collectives.hpp>
#include <boost/mpi/communicator.hpp>

#include <algorithm>
#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace ScriptInterface {
namespace Particles {

#ifdef EXCLUSIONS
/**
 * @brief Use the bond topology to automatically add exclusions between
 * particles that are up to @c n_bonds_max bonds apart in a chain.
 */
static void auto_exclusions(boost::mpi::communicator const &comm,
                            int const n_bonds_max) {
  // bookkeeping of particle exclusions, with their n-bond distance
  std::unordered_map<int, std::vector<std::pair<int, int>>> partners;
  std::vector<int> bonded_pairs;

  auto &system = ::System::get_system();
  auto &cell_structure = *system.cell_structure;

  // determine initial connectivity
  for (auto const &p : cell_structure.local_particles()) {
    auto const pid1 = p.id();
    for (auto const bond : p.bonds()) {
      if (bond.partner_ids().size() == 1u) {
        auto const pid2 = bond.partner_ids()[0];
        if (pid1 != pid2) {
          bonded_pairs.emplace_back(pid1);
          bonded_pairs.emplace_back(pid2);
        }
      }
    }
  }

  Utils::Mpi::gather_buffer(bonded_pairs, comm);

  if (comm.rank() == 0) {
    auto const add_partner = [&partners](int pid1, int pid2, int n_bonds) {
      if (pid2 == pid1)
        return;
      for (auto const &partner : partners[pid1])
        if (partner.first == pid2)
          return;
      partners[pid1].emplace_back(pid2, n_bonds);
    };

    for (auto it = bonded_pairs.begin(); it != bonded_pairs.end(); it += 2) {
      add_partner(it[0], it[1], 1);
      add_partner(it[1], it[0], 1);
    }

    // determine transient connectivity
    for (int iteration = 1; iteration < n_bonds_max; iteration++) {
      std::vector<int> pids;
      for (auto const &kv : partners) {
        pids.emplace_back(kv.first);
      }
      for (auto const pid1 : pids) {
        // loop over partners (counter-based loops due to iterator invalidation)
        // NOLINTNEXTLINE(modernize-loop-convert)
        for (std::size_t i = 0u; i < partners[pid1].size(); ++i) {
          auto const [pid2, dist21] = partners[pid1][i];
          if (dist21 > n_bonds_max)
            continue;
          // loop over all partners of the partner
          // NOLINTNEXTLINE(modernize-loop-convert)
          for (std::size_t j = 0u; j < partners[pid2].size(); ++j) {
            auto const [pid3, dist32] = partners[pid2][j];
            auto const dist31 = dist32 + dist21;
            if (dist31 > n_bonds_max)
              continue;
            add_partner(pid1, pid3, dist31);
            add_partner(pid3, pid1, dist31);
          }
        }
      }
    }
  }

  boost::mpi::broadcast(comm, partners, 0);
  for (auto const &kv : partners) {
    auto const pid1 = kv.first;
    auto const &partner_list = kv.second;
    for (auto const &partner : partner_list) {
      auto const pid2 = partner.first;
      if (auto p1 = cell_structure.get_local_particle(pid1)) {
        add_exclusion(*p1, pid2);
      }
      if (auto p2 = cell_structure.get_local_particle(pid2)) {
        add_exclusion(*p2, pid1);
      }
    }
  }
  system.on_particle_change();
}
#endif // EXCLUSIONS

Variant ParticleList::do_call_method(std::string const &name,
                                     VariantMap const &params) {
#ifdef EXCLUSIONS
  if (name == "auto_exclusions") {
    auto const distance = get_value<int>(params, "distance");
    auto_exclusions(context()->get_comm(), distance);
    return {};
  }
#endif // EXCLUSIONS
  if (name == "get_highest_particle_id") {
    return get_maximal_particle_id();
  }
  if (name == "clear") {
    remove_all_particles();
    return {};
  }
  if (not context()->is_head_node()) {
    return {};
  }
  if (name == "by_id") {
    return std::dynamic_pointer_cast<ParticleHandle>(
        context()->make_shared("Particles::ParticleHandle",
                               {{"id", get_value<int>(params, "p_id")},
                                {"__cell_structure", m_cell_structure.lock()},
                                {"__bonded_ias", m_bonded_ias.lock()}}));
  }
  if (name == "by_ids") {
    return context()->make_shared(
        "Particles::ParticleSlice",
        {{"id_selection", get_value<std::vector<int>>(params, "id_selection")},
         {"__cell_structure", m_cell_structure.lock()},
         {"__bonded_ias", m_bonded_ias.lock()}});
  }
  if (name == "get_n_part") {
    return get_n_part();
  }
  if (name == "get_particle_ids") {
    return get_particle_ids();
  }
  if (name == "particle_exists") {
    return particle_exists(get_value<int>(params, "p_id"));
  }
  if (name == "add_particle") {
    assert(not params.contains("bonds"));
    VariantMap local_params = params;
    local_params["__cell_structure"] = m_cell_structure.lock();
    local_params["__bonded_ias"] = m_bonded_ias.lock();
    auto so = std::dynamic_pointer_cast<ParticleHandle>(
        context()->make_shared("Particles::ParticleHandle", local_params));
#ifdef EXCLUSIONS
    if (params.count("exclusions")) {
      so->call_method("set_exclusions", {{"p_ids", params.at("exclusions")}});
    }
#endif // EXCLUSIONS
    return so->get_parameter("id");
  }
  return {};
}

} // namespace Particles
} // namespace ScriptInterface
