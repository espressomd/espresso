/*
 * Copyright (C) 2021-2022 The ESPResSo project
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

#ifndef ESPRESSO_SRC_SCRIPT_INTERFACE_WALBERLA_LATTICE_WALBERLA_HPP
#define ESPRESSO_SRC_SCRIPT_INTERFACE_WALBERLA_LATTICE_WALBERLA_HPP

#include "config/config.hpp"

#ifdef WALBERLA

#include <walberla_bridge/LatticeWalberla.hpp>
#include <walberla_bridge/lattice_boltzmann/lb_walberla_init.hpp>

#include "core/errorhandling.hpp"
#include "core/grid.hpp"

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"

#include <utils/Vector.hpp>

#include <cassert>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

namespace ScriptInterface::walberla {

class LatticeWalberla : public AutoParameters<LatticeWalberla> {
  std::shared_ptr<::LatticeWalberla> m_lattice;
  int m_n_ghost_layers;
  double m_agrid;

public:
  LatticeWalberla() {
    add_parameters(
        {{"agrid", AutoParameter::read_only, [this]() { return m_agrid; }},
         {"n_ghost_layers", AutoParameter::read_only,
          [this]() { return m_n_ghost_layers; }}});
  }

  void do_construct(VariantMap const &args) override {
    auto const box_size = box_geo.length();
    auto const agrid = get_value<double>(args, "agrid");
    auto const n_ghost_layers = get_value<int>(args, "n_ghost_layers");
    assert(agrid > 0.);
    assert(n_ghost_layers >= 0);
    m_agrid = agrid;
    m_n_ghost_layers = n_ghost_layers;

    try {
      auto const grid_dimensions = ::calc_grid_dimensions(box_size, agrid);
      m_lattice = std::make_shared<::LatticeWalberla>(
          grid_dimensions, node_grid, static_cast<unsigned>(n_ghost_layers));
    } catch (const std::exception &e) {
      runtimeErrorMsg() << "LatticeWalberla failed: " << e.what();
    }
  }

  std::shared_ptr<::LatticeWalberla> lattice() { return m_lattice; }
  std::shared_ptr<const ::LatticeWalberla> lattice() const { return m_lattice; }
};

} // namespace ScriptInterface::walberla

#endif // WALBERLA
#endif
