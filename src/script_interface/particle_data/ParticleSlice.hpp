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

#include "ParticleHandle.hpp"

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/Variant.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"
#include "script_interface/cell_system/CellSystem.hpp"
#include "script_interface/get_value.hpp"
#include "script_interface/interactions/BondedInteractions.hpp"

#include "core/system/System.hpp"

#include <utils/mpi/gather_buffer.hpp>

#include <algorithm>
#include <cassert>
#include <functional>
#include <memory>
#include <ranges>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace ScriptInterface {
namespace Particles {

struct SetParticleParametersVisitor {
  void operator()(std::vector<int> const &, std::string const &,
                  auto const &values, Context *,
                  std::shared_ptr<CellSystem::CellSystem>,
                  std::shared_ptr<Interactions::BondedInteractions>) const {
    throw Exception("Values must be of type vector, got " +
                    detail::demangle::simplify_symbol(&values));
  }
  template <typename T>
  void operator()(
      std::vector<int> const &pids, std::string const &param_name,
      std::vector<T> const &values, Context *context,
      std::shared_ptr<CellSystem::CellSystem> cell_structure,
      std::shared_ptr<Interactions::BondedInteractions> bonded_ias) const {
    auto so = std::dynamic_pointer_cast<ParticleModifier>(context->make_shared(
        "Particles::ParticleModifier", {{"id", -1},
                                        {"__cell_structure", cell_structure},
                                        {"__bonded_ias", bonded_ias}}));
    for (std::size_t i = 0; i < pids.size(); ++i) {
      so->set_pid(pids[i]);
      so->do_set_parameter(param_name, values[i]);
    }
  }
};

template <typename T>
inline auto
get_particles_properties(std::vector<int> const &pids,
                         std::function<T(Particle const &)> const &getter,
                         Context *context,
                         CellStructure const &cell_structure) {

  using value_type =
      std::conditional_t<Variant::has_type<std::vector<T>>::value, T, Variant>;
  std::vector<value_type> result;

  auto const n_ranks = static_cast<std::size_t>(context->get_comm().size());
  auto const size_hint = pids.size() / n_ranks;
  std::vector<std::pair<int, T>> parameters;
  parameters.reserve(size_hint);
  for (auto const &pid : pids) {
    auto const p = cell_structure.get_local_particle(pid);
    if (p and not p->is_ghost()) {
      parameters.emplace_back(pid, getter(*p));
    }
  }

  // collect values from all nodes
  Utils::Mpi::gather_buffer(parameters, context->get_comm(), 0);
  if (!context->is_head_node()) {
    return result;
  }

  // sort values by particle id to retain original order
  auto const projector = [](auto const &pair) { return pair.first; };
  std::ranges::sort(parameters, std::less<int>{}, projector);
  assert(std::ranges::equal(pids, parameters | std::views::keys) &&
         "Missing or duplicate particle ids");

  result.reserve(pids.size());
  for (auto const &value : parameters | std::views::values) {
    result.emplace_back(std::move(value));
  }
  return result;
}

class ParticleSlice : public AutoParameters<ParticleSlice> {
  std::vector<int> m_id_selection;
  int m_chunk_size;
  std::weak_ptr<CellSystem::CellSystem> m_cell_structure;
  std::weak_ptr<Interactions::BondedInteractions> m_bonded_ias;
  std::weak_ptr<::System::System> m_system;
  /** @brief Data structure to store names of parameters with special setters.
   */
  std::set<std::string_view> const m_special_parameters{
      "pos",        "type", "bonds",
#ifdef ESPRESSO_ELECTROSTATICS
      "q",
#endif // ESPRESSO_ELECTROSTATICS
#ifdef ESPRESSO_EXCLUSIONS
      "exclusions",
#endif // ESPRESSO_EXCLUSIONS
  };

  auto get_cell_structure() const {
    auto cell_structure_ptr = m_cell_structure.lock();
    assert(cell_structure_ptr != nullptr);
    auto &cell_structure = cell_structure_ptr->get_cell_structure();
    return &cell_structure;
  }

  auto get_system() const {
    auto ptr = m_system.lock();
    assert(ptr != nullptr);
    return ptr;
  }

public:
  ParticleSlice() {
    add_parameters({
        {"chunk_size", AutoParameter::read_only,
         [this]() { return m_chunk_size; }},
        {"id_selection", AutoParameter::read_only,
         [this]() { return m_id_selection; }},
    });
  }

  void do_construct(VariantMap const &params) override;

  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override;

  void attach(std::weak_ptr<::System::System> system) {
    assert(m_system.expired());
    m_system = system;
  }
};

} // namespace Particles
} // namespace ScriptInterface
