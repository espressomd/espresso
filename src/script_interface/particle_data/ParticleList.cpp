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

#include "script_interface/ObjectState.hpp"
#include "script_interface/ScriptInterface.hpp"

#include "core/cells.hpp"
#include "core/event.hpp"
#include "core/exclusions.hpp"
#include "core/particle_node.hpp"

#include <utils/Vector.hpp>
#include <utils/mpi/gather_buffer.hpp>
#include <utils/serialization/pack.hpp>

#include <boost/mpi/collectives.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/range/algorithm.hpp>

#include <memory>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace ScriptInterface {
namespace Particles {

#ifdef EXCLUSIONS
static void set_exclusions(ParticleHandle &p, Variant const &exclusions) {
  p.call_method("set_exclusions", {{"p_ids", exclusions}});
}
#endif // EXCLUSIONS

static void set_bonds(ParticleHandle &p, Variant const &bonds) {
  auto const bond_list_flat = get_value<std::vector<std::vector<int>>>(bonds);
  for (auto const &bond_flat : bond_list_flat) {
    auto const bond_id = bond_flat[0];
    auto const part_id =
        std::vector<int>{bond_flat.begin() + 1, bond_flat.end()};
    p.call_method("add_bond",
                  {{"bond_id", bond_id}, {"part_id", std::move(part_id)}});
  }
}

std::string ParticleList::get_internal_state() const {
  auto const p_ids = get_particle_ids();
  std::vector<std::string> object_states(p_ids.size());

  boost::transform(p_ids, object_states.begin(), [this](auto const p_id) {
    auto p_obj =
        context()->make_shared("Particles::ParticleHandle", {{"id", p_id}});
    auto &p_handle = dynamic_cast<ParticleHandle &>(*p_obj);
    auto const packed_state = p_handle.serialize();
    // custom particle serialization
    auto state = Utils::unpack<ObjectState>(packed_state);
    state.name = "Particles::ParticleHandle";
    auto const bonds_view = p_handle.call_method("get_bonds_view", {});
    state.params.emplace_back(std::string{"bonds"}, pack(bonds_view));
#ifdef EXCLUSIONS
    auto const exclusions = p_handle.call_method("get_exclusions", {});
    state.params.emplace_back(std::string{"exclusions"}, pack(exclusions));
#endif // EXCLUSIONS
    state.params.emplace_back(std::string{"__cpt_sentinel"}, pack(None{}));
    return Utils::pack(state);
  });

  return Utils::pack(object_states);
}

void ParticleList::set_internal_state(std::string const &state) {
  auto const object_states = Utils::unpack<std::vector<std::string>>(state);
#ifdef EXCLUSIONS
  std::unordered_map<int, Variant> exclusions = {};
#endif // EXCLUSIONS
  std::unordered_map<int, Variant> bonds = {};

  for (auto const &packed_object : object_states) {
    auto state = Utils::unpack<ObjectState>(packed_object);
    VariantMap params = {};
    for (auto const &kv : state.params) {
      params[kv.first] = unpack(kv.second, {});
    }
    auto const p_id = get_value<int>(params.at("id"));
    bonds[p_id] = params.extract("bonds").mapped();
#ifdef EXCLUSIONS
    exclusions[p_id] = params.extract("exclusions").mapped();
#endif // EXCLUSIONS
    context()->make_shared("Particles::ParticleHandle", params);
  }

  for (auto const p_id : get_particle_ids()) {
    auto p_obj =
        context()->make_shared("Particles::ParticleHandle", {{"id", p_id}});
    auto &p_handle = dynamic_cast<ParticleHandle &>(*p_obj);
    set_bonds(p_handle, bonds[p_id]);
#ifdef EXCLUSIONS
    set_exclusions(p_handle, exclusions[p_id]);
#endif // EXCLUSIONS
  }
}

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

  // determine initial connectivity
  for (auto const &p : ::cell_structure.local_particles()) {
    auto const pid1 = p.id();
    for (auto const bond : p.bonds()) {
      if (bond.partner_ids().size() == 1) {
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
        for (int i = 0; i < partners[pid1].size(); ++i) {
          auto const [pid2, dist21] = partners[pid1][i];
          if (dist21 > n_bonds_max)
            continue;
          // loop over all partners of the partner
          // NOLINTNEXTLINE(modernize-loop-convert)
          for (int j = 0; j < partners[pid2].size(); ++j) {
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
      if (auto p1 = ::cell_structure.get_local_particle(pid1)) {
        add_exclusion(*p1, pid2);
      }
      if (auto p2 = ::cell_structure.get_local_particle(pid2)) {
        add_exclusion(*p2, pid1);
      }
    }
  }
  on_particle_change();
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
    assert(params.count("bonds") == 0);
    auto obj = context()->make_shared("Particles::ParticleHandle", params);
    auto &p_handle = dynamic_cast<ParticleHandle &>(*obj);
#ifdef EXCLUSIONS
    if (params.count("exclusions")) {
      set_exclusions(p_handle, params.at("exclusions"));
    }
#endif // EXCLUSIONS
    return p_handle.get_parameter("id");
  }
  return {};
}

} // namespace Particles
} // namespace ScriptInterface
