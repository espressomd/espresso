/*
 * Copyright (C) 2021-2023 The ESPResSo project
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

#include "config/config.hpp"

#ifdef WALBERLA

#include "core/grid.hpp"

#include <walberla_bridge/LatticeWalberla.hpp>

#include <script_interface/ScriptInterface.hpp>
#include <script_interface/auto_parameters/AutoParameters.hpp>

#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>

namespace ScriptInterface::walberla {

class LatticeWalberla : public AutoParameters<LatticeWalberla> {
  std::shared_ptr<::LatticeWalberla> m_lattice;
  double m_agrid;

public:
  LatticeWalberla() {
    add_parameters({
        {"agrid", AutoParameter::read_only, [this]() { return m_agrid; }},
        {"n_ghost_layers", AutoParameter::read_only,
         [this]() { return static_cast<int>(m_lattice->get_ghost_layers()); }},
        {"shape", AutoParameter::read_only,
         [this]() { return m_lattice->get_grid_dimensions(); }},
    });
  }

  void do_construct(VariantMap const &args) override {
    m_agrid = get_value<double>(args, "agrid");
    auto const n_ghost_layers = get_value<int>(args, "n_ghost_layers");

    context()->parallel_try_catch([&]() {
      if (m_agrid <= 0.) {
        throw std::domain_error("Parameter 'agrid' must be > 0");
      }
      if (n_ghost_layers < 0) {
        throw std::domain_error("Parameter 'n_ghost_layers' must be >= 0");
      }
      auto const box_size = ::box_geo.length();
      auto const grid_dim =
          ::LatticeWalberla::calc_grid_dimensions(box_size, m_agrid);
      m_lattice = std::make_shared<::LatticeWalberla>(
          grid_dim, node_grid, static_cast<unsigned int>(n_ghost_layers));
    });
  }

  std::shared_ptr<::LatticeWalberla> lattice() { return m_lattice; }
  std::shared_ptr<const ::LatticeWalberla> lattice() const { return m_lattice; }
};

} // namespace ScriptInterface::walberla

#endif // WALBERLA
