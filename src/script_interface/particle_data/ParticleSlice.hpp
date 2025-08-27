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
#include "script_interface/auto_parameters/AutoParameters.hpp"
#include "script_interface/cell_system/CellSystem.hpp"
#include "script_interface/interactions/BondedInteractions.hpp"

#include "core/system/System.hpp"

#include <memory>
#include <string>
#include <vector>

namespace ScriptInterface {
namespace Particles {

void set_particles_bonds(
    std::vector<int> pids, std::vector<std::vector<int>> all_bonds_ids,
    std::vector<std::vector<std::vector<int>>> all_bonds_partner_ids,
    ::CellStructure &cell_structure, Context *context,
    std::shared_ptr<::System::System> system);

#ifdef ESPRESSO_EXCLUSIONS
void set_particles_exclusions(std::vector<int> pids,
                              std::vector<std::vector<int>> exclusion_lists,
                              ::CellStructure &cell_structure, Context *context,
                              std::shared_ptr<::System::System> system);
#endif // ESPRESSO_EXCLUSIONS

void set_particles_positions(std::vector<int> pids,
                             std::vector<Utils::Vector3d> positions);

void set_particles_types(std::vector<int> pids, std::vector<int> types,
                         boost::mpi::communicator &comm,
                         CellStructure &cell_structure,
                         std::shared_ptr<::System::System> system);

#ifdef ESPRESSO_ELECTROSTATICS
void set_particles_charges(
    std::vector<int> pids,
    std::variant<std::vector<int>, std::vector<double>> const &charges,
    boost::mpi::communicator const &comm, CellStructure &cell_structure,
    std::shared_ptr<::System::System> system);
#endif // ESPRESSO_ELECTROSTATICS

struct SetParticleParametersVisitor {
  template <typename T>
  void operator()(
      std::vector<int> const pids, std::basic_string<char> const param_name,
      std::vector<T> values, Context *context,
      std::shared_ptr<CellSystem::CellSystem> cell_structure,
      std::shared_ptr<Interactions::BondedInteractions> bonded_ias) const {
    for (std::size_t i = 0; i < pids.size(); ++i) {
      auto const pid = pids[i];
      context
          ->make_shared("Particles::ParticleHandle",
                        {{"id", pid},
                         {"__cell_structure", cell_structure},
                         {"__bonded_ias", bonded_ias}})
          ->do_set_parameter(param_name, values[i]);
    } // Particle id loop
  };
  void operator()(
      std::vector<int> const pids, std::basic_string<char> const param_name,
      auto values, Context *context,
      std::shared_ptr<CellSystem::CellSystem> cell_structure,
      std::shared_ptr<Interactions::BondedInteractions> bonded_ias) const {
    throw Exception("Values must be of type vector.");
  }
};

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
