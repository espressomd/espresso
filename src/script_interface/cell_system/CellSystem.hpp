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

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"
#include "script_interface/system/Leaf.hpp"
#include "script_interface/system/System.hpp"

#include "core/cell_system/CellStructure.hpp"
#include "core/cell_system/CellStructureType.hpp"
#include "core/cell_system/HybridDecomposition.hpp"
#include "core/cell_system/RegularDecomposition.hpp"

#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace ScriptInterface {
namespace Particles {
class ParticleHandle;
class ParticleSlice;
} // namespace Particles
namespace CellSystem {

class CellSystem : public AutoParameters<CellSystem, System::Leaf> {
  std::unordered_map<CellStructureType, std::string> const cs_type_to_name = {
      {CellStructureType::REGULAR, "regular_decomposition"},
      {CellStructureType::NSQUARE, "n_square"},
      {CellStructureType::HYBRID, "hybrid_decomposition"},
  };

  std::unordered_map<std::string, CellStructureType> const cs_name_to_type = {
      {"regular_decomposition", CellStructureType::REGULAR},
      {"n_square", CellStructureType::NSQUARE},
      {"hybrid_decomposition", CellStructureType::HYBRID},
  };

  std::shared_ptr<::CellStructure> m_cell_structure;
  std::unique_ptr<VariantMap> m_params;

  void on_bind_system(::System::System &system) override {
    m_cell_structure = system.cell_structure;
    m_cell_structure->bind_system(m_system.lock());
    auto const &params = *m_params;
    if (not params.empty()) {
      auto const cs_name = get_value<std::string>(params, "decomposition_type");
      auto const cs_type = cs_name_to_type.at(cs_name);
      initialize(cs_type, params);
      do_set_parameter("skin", params.at("skin"));
      do_set_parameter("node_grid", params.at("node_grid"));
    }
    m_params.reset();
  }

public:
  CellSystem();

  void do_construct(VariantMap const &params) override {
    m_params = std::make_unique<VariantMap>(params);
  }

  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override;

  auto &get_cell_structure() const { return *m_cell_structure; }

  void configure(Particles::ParticleHandle &);
  void configure(Particles::ParticleSlice &);

private:
  /**
   * @brief Resort the particles.
   *
   * This function resorts the particles on the nodes.
   *
   * @param global_flag   If true a global resort is done, if false particles
   *                      are only exchanges between neighbors.
   * @return The number of particles on the nodes after the resort.
   */
  std::vector<int> mpi_resort_particles(bool global_flag) const;

  void initialize(CellStructureType const &cs_type, VariantMap const &params);

  auto const &get_regular_decomposition() const {
    return dynamic_cast<RegularDecomposition const &>(
        std::as_const(get_cell_structure()).decomposition());
  }

  auto const &get_hybrid_decomposition() const {
    return dynamic_cast<HybridDecomposition const &>(
        std::as_const(get_cell_structure()).decomposition());
  }
};

} // namespace CellSystem
} // namespace ScriptInterface
