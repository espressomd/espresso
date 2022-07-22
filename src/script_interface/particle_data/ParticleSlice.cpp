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

#include "script_interface/ScriptInterface.hpp"

#include "core/particle_node.hpp"

#include <utils/Span.hpp>

#include <string>
#include <vector>

namespace ScriptInterface {
namespace Particles {

void ParticleSlice::do_construct(VariantMap const &params) {
  m_id_selection = get_value<std::vector<int>>(params, "id_selection");
  m_chunk_size = get_value_or<int>(params, "prefetch_chunk_size", 10000);
}

Variant ParticleSlice::do_call_method(std::string const &name,
                                      VariantMap const &params) {
  if (name == "prefetch_particle_data") {
    auto p_ids = get_value<std::vector<int>>(params, "chunk");
    prefetch_particle_data(Utils::Span<int>(p_ids));
  } else if (name == "particle_exists") {
    return particle_exists(get_value<int>(params, "p_id"));
  }
  return {};
}

} // namespace Particles
} // namespace ScriptInterface
