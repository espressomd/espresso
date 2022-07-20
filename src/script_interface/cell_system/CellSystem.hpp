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

#ifndef ESPRESSO_SRC_SCRIPT_INTERFACE_CELL_SYSTEM_CELL_SYSTEM_HPP
#define ESPRESSO_SRC_SCRIPT_INTERFACE_CELL_SYSTEM_CELL_SYSTEM_HPP

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"

#include "cell_system/CellStructureType.hpp"

#include <string>
#include <unordered_map>

namespace ScriptInterface {
namespace CellSystem {

class CellSystem : public AutoParameters<CellSystem> {
  std::unordered_map<CellStructureType, std::string> const cs_type_to_name = {
      {CellStructureType::CELL_STRUCTURE_REGULAR, "regular_decomposition"},
      {CellStructureType::CELL_STRUCTURE_NSQUARE, "n_square"},
      {CellStructureType::CELL_STRUCTURE_HYBRID, "hybrid_decomposition"},
  };

  std::unordered_map<std::string, CellStructureType> const cs_name_to_type = {
      {"regular_decomposition", CellStructureType::CELL_STRUCTURE_REGULAR},
      {"n_square", CellStructureType::CELL_STRUCTURE_NSQUARE},
      {"hybrid_decomposition", CellStructureType::CELL_STRUCTURE_HYBRID},
  };

public:
  CellSystem();

  void do_construct(VariantMap const &params) override {
    // the core cell structure has a default constructor; params is empty
    // during default construction for the python class, and not empty
    // during checkpointing
    if (params.count("decomposition_type")) {
      auto const cs_name = get_value<std::string>(params, "decomposition_type");
      auto const cs_type = cs_name_to_type.at(cs_name);
      initialize(cs_type, params);
      do_set_parameter("skin", params.at("skin"));
      do_set_parameter("node_grid", params.at("node_grid"));
    }
  }

  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override;

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

  void initialize(CellStructureType const &cs_type,
                  VariantMap const &params) const;
};

} // namespace CellSystem
} // namespace ScriptInterface

#endif
