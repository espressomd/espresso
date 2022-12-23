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

#include "core/particle_node.hpp"

#include <utils/Vector.hpp>
#include <utils/serialization/pack.hpp>

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
  p.do_call_method("set_exclusions", {{"p_ids", exclusions}});
}
#endif // EXCLUSIONS

static void set_bonds(ParticleHandle &p, Variant const &bonds) {
  auto const bond_list_flat = get_value<std::vector<std::vector<int>>>(bonds);
  for (auto const &bond_flat : bond_list_flat) {
    auto const bond_id = bond_flat[0];
    auto const part_id =
        std::vector<int>{bond_flat.begin() + 1, bond_flat.end()};
    p.do_call_method("add_bond",
                     {{"bond_id", bond_id}, {"part_id", std::move(part_id)}});
  }
}

std::string ParticleList::get_internal_state() const {
  auto const p_ids = get_particle_ids();
  std::vector<std::string> object_states(p_ids.size());

  boost::transform(p_ids, object_states.begin(), [](auto const p_id) {
    ParticleHandle p_handle{};
    p_handle.do_construct({{"id", p_id}});
    auto const packed_state = p_handle.serialize();
    // custom particle serialization
    auto state = Utils::unpack<ObjectState>(packed_state);
    state.name = "Particles::ParticleHandle";
    auto const bonds_view = p_handle.call_method("get_bonds_view", {});
    state.params.emplace_back(
        std::pair<std::string, PackedVariant>{"bonds", pack(bonds_view)});
#ifdef EXCLUSIONS
    auto const exclusions = p_handle.call_method("get_exclusions", {});
    state.params.emplace_back(
        std::pair<std::string, PackedVariant>{"exclusions", pack(exclusions)});
#endif // EXCLUSIONS
    state.params.emplace_back(
        std::pair<std::string, PackedVariant>{"__cpt_sentinel", pack(None{})});
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
    auto const state = Utils::unpack<ObjectState>(packed_object);
    auto o = std::dynamic_pointer_cast<ParticleHandle>(
        ObjectHandle::deserialize(packed_object, *ObjectHandle::context()));
    auto const p_id = get_value<int>(o->get_parameter("id"));
    for (auto const &kv : state.params) {
      if (kv.first == "bonds") {
        bonds[p_id] = unpack(kv.second, {});
      }
#ifdef EXCLUSIONS
      else if (kv.first == "exclusions") {
        exclusions[p_id] = unpack(kv.second, {});
      }
#endif // EXCLUSIONS
    }
  }

  for (auto const p_id : get_particle_ids()) {
    ParticleHandle p_handle{};
    p_handle.do_construct({{"id", p_id}});
    set_bonds(p_handle, bonds[p_id]);
#ifdef EXCLUSIONS
    set_exclusions(p_handle, exclusions[p_id]);
#endif // EXCLUSIONS
  }
}

Variant ParticleList::do_call_method(std::string const &name,
                                     VariantMap const &params) {
  if (name == "get_n_part") {
    return get_n_part();
  }
  if (name == "get_particle_ids") {
    return get_particle_ids();
  }
  if (name == "particle_exists") {
    return particle_exists(get_value<int>(params, "p_id"));
  }
  if (name == "clear") {
    remove_all_particles();
  } else if (name == "add_particle") {
    assert(params.count("bonds") == 0);
    // sanitize particle properties
#ifdef DIPOLES
    if (params.count("dip") and params.count("dipm")) {
      throw std::invalid_argument("Contradicting attributes: 'dip' and 'dipm'. \
Setting 'dip' is sufficient as the length of the vector defines the scalar \
dipole moment.");
    }
    if (params.count("dip") and params.count("quat")) {
      throw std::invalid_argument("Contradicting attributes: 'dip' and 'quat'. \
Setting 'dip' overwrites the rotation of the particle around the dipole axis. \
Set attribute 'quat' together with 'dipm' (scalar dipole moment) instead.");
    }
#endif // DIPOLES
    ParticleHandle p_handle{};
    p_handle.do_construct(params);
#ifdef EXCLUSIONS
    if (params.count("exclusions")) {
      set_exclusions(p_handle, params.at("exclusions"));
    }
#endif // EXCLUSIONS
    return p_handle.get_parameter("id");
  }
  if (name == "get_highest_particle_id") {
    return get_maximal_particle_id();
  }
  return {};
}

} // namespace Particles
} // namespace ScriptInterface
