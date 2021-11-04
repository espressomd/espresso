/*
 * Copyright (C) 2021 The ESPResSo project
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
#ifndef SCRIPT_INTERFACE_WALBERLA_FLUIDWALBERLA_HPP
#define SCRIPT_INTERFACE_WALBERLA_FLUIDWALBERLA_HPP

#include "LatticeWalberla.hpp"

#include "core/grid_based_algorithms/lb_walberla_instance.hpp"
#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"

#include <utils/Vector.hpp>

#include <memory>
#include <string>

namespace ScriptInterface::walberla {

class FluidWalberla : public AutoParameters<FluidWalberla> {
public:
  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override {
    if (name == "init_lb_walberla") {
      auto const lb_visc = get_value<double>(params, "visc");
      auto const lb_dens = get_value<double>(params, "dens");
      auto const lb_kT = get_value<double>(params, "kT");
      auto const lb_seed = get_value<int>(params, "seed");
      auto const tau = get_value<double>(params, "tau");
      auto const agrid = get_value<double>(params, "agrid");
      auto const lattice =
          get_value<std::shared_ptr<LatticeWalberla>>(params, "lattice")
              ->lattice();
      mpi_init_lb_walberla_local(lattice, lb_visc, lb_dens, agrid, tau, lb_kT,
                                 lb_seed);
      return {};
    }

    return {};
  }
};

} // namespace ScriptInterface::walberla

#endif
