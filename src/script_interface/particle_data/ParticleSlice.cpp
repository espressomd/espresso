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

#include "ParticleSlice.hpp"
#include "ParticleHandle.hpp"

#include "script_interface/ScriptInterface.hpp"

#include "core/particle_node.hpp"

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace ScriptInterface {
namespace Particles {

void ParticleSlice::do_construct(VariantMap const &params) {
  if (params.contains("__cell_structure")) {
    auto so = get_value<std::shared_ptr<CellSystem::CellSystem>>(
        params, "__cell_structure");
    so->configure(*this);
    m_cell_structure = so;
  }
  if (params.contains("__bonded_ias")) {
    m_bonded_ias = get_value<std::shared_ptr<Interactions::BondedInteractions>>(
        params, "__bonded_ias");
  }
  m_id_selection = get_value<std::vector<int>>(params, "id_selection");
  m_chunk_size = get_value_or<int>(params, "prefetch_chunk_size", 10000);
  if (not context()->is_head_node()) {
    return;
  }
  for (auto const pid : m_id_selection) {
    if (not particle_exists(pid)) {
      throw std::out_of_range("Particle does not exist: " +
                              std::to_string(pid));
    }
  }
}

Variant ParticleSlice::do_call_method(std::string const &name,
                                      VariantMap const &params) {
  if (not context()->is_head_node()) {
    return {};
  }
  if (name == "prefetch_particle_data") {
    auto p_ids = get_value<std::vector<int>>(params, "chunk");
    prefetch_particle_data(p_ids);
    return {};
  }
  if (name == "get_particle") {
    return context()->make_shared(
        "Particles::ParticleHandle",
        {{"id", get_value<int>(params, "p_id")},
         {"__cell_structure", m_cell_structure.lock()},
         {"__bonded_ias", m_bonded_ias.lock()}});
  }
  return {};
}

} // namespace Particles
} // namespace ScriptInterface
