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
#ifndef SCRIPT_INTERFACE_WALBERLA_FLUID_ROUTINES_WALBERLA_HPP
#define SCRIPT_INTERFACE_WALBERLA_FLUID_ROUTINES_WALBERLA_HPP

#include "config.hpp"

#ifdef LB_WALBERLA

#include <walberla_bridge/LBWalberlaBase.hpp>

#include "core/errorhandling.hpp"
#include "core/grid_based_algorithms/lb_interface.hpp"
#include "core/grid_based_algorithms/lb_walberla_instance.hpp"
#include "core/grid_based_algorithms/lb_walberla_interface.hpp"
#include "core/integrate.hpp"

#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"

#include <utils/Vector.hpp>
#include <utils/constants.hpp>
#include <utils/math/int_pow.hpp>

#include <memory>
#include <stdexcept>
#include <string>

namespace ScriptInterface::walberla {

class FluidNodeWalberla : public AutoParameters<FluidNodeWalberla> {
  // TODO WALBERLA: revert commit 2f0c490b8e1bb4ab3e to use a weak_ptr
  std::shared_ptr<::LBWalberlaBase> m_lb_fluid;
  Utils::Vector3i m_index;
  Utils::Vector3i m_grid_size;
  double m_conv_dens;
  double m_conv_press;
  double m_conv_force;
  double m_conv_velocity;

public:
  FluidNodeWalberla();

  void do_construct(VariantMap const &params) override {
    try {
      m_lb_fluid = ::lb_walberla();
      auto const &lb_params = ::lb_walberla_params();
      auto const tau = lb_params->get_tau();
      auto const agrid = lb_params->get_agrid();
      m_grid_size = m_lb_fluid->lattice().get_grid_dimensions();
      m_index = get_value<Utils::Vector3i>(params, "index");
      if (not(is_index_valid(m_index))) {
        throw std::out_of_range("Index error");
      }
      m_conv_dens = Utils::int_pow<3>(agrid);
      m_conv_press = Utils::int_pow<1>(agrid) * Utils::int_pow<2>(tau);
      m_conv_force = Utils::int_pow<2>(tau) / Utils::int_pow<1>(agrid);
      m_conv_velocity = Utils::int_pow<1>(tau) / Utils::int_pow<1>(agrid);
    } catch (const std::exception &e) {
      runtimeErrorMsg() << "LatticeWalberla failed: " << e.what();
      m_lb_fluid.reset();
    }
  }

  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override {
    if (name == "override_index") {
      // this hidden feature is used to iterate a LB slice without
      // rebuilding a FluidNodeWalberla for each node in the slice
      auto const index = get_value<Utils::Vector3i>(params, "index");
      if (not(is_index_valid(index))) {
        return ES_ERROR;
      }
      m_index = index;
      return ES_OK;
    }

    return {};
  }

private:
  inline bool is_index_valid(Utils::Vector3i const &index) const {
    return index < m_grid_size and index >= Utils::Vector3i{};
  }
};
} // namespace ScriptInterface::walberla

#endif // LB_WALBERLA
#endif
