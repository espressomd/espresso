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
#include "utils/mpi/gather_buffer.hpp"

#include <functional>
#include <memory>
#include <string>
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
inline std::vector<Variant>
get_particles_properties(std::vector<int> const &pids,
                         std::function<T(Particle const &)> const &getter,
                         Context *context, CellStructure &cell_structure,
                         ::System::System const &system) {
  std::vector<Variant> result;
  result.reserve(pids.size());
  std::vector<std::pair<int, T>> index_value_vec;
  for (auto const &pid : pids) {
    auto const p = cell_structure.get_local_particle(pid);
    if (p and not p->is_ghost()) {
      index_value_vec.emplace_back(pid, getter(*p));
    }
  }
  // Collect values from all nodes
  Utils::Mpi::gather_buffer(index_value_vec, context->get_comm(), 0);
  if (!context->is_head_node()) {
    return {};
  }

  // Sort results by particle ids to retain original order
  std::ranges::sort(index_value_vec, [](const auto &pair1, const auto &pair2) {
    return pair1.first < pair2.first;
  });
  assert(std::ranges::equal(
             pids, index_value_vec,
             [](int const &pid1, std::pair<int, T> const &index_value_pair) {
               return pid1 == index_value_pair.first;
             }) &&
         "Unexpected returned values.");
  for (auto const &index_value_pair : index_value_vec) {
    result.emplace_back(std::move(index_value_pair.second));
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
  std::set<std::string> m_special_parameters{
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
