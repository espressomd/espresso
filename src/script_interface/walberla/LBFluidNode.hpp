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

#include "LBFluid.hpp"

#include "LatticeIndices.hpp"

#include <script_interface/ScriptInterface.hpp>
#include <script_interface/auto_parameters/AutoParameters.hpp>

#include <walberla_bridge/lattice_boltzmann/LBWalberlaBase.hpp>

#include <utils/Vector.hpp>
#include <utils/math/int_pow.hpp>

#include <cassert>
#include <memory>
#include <stdexcept>
#include <string>

namespace ScriptInterface::walberla {

class LBFluidNode : public AutoParameters<LBFluidNode, LatticeIndices> {
  std::shared_ptr<::LBWalberlaBase> m_lb_fluid;
  Utils::Vector3i m_index;
  Utils::Vector3i m_grid_size;
  double m_conv_dens;
  double m_conv_press;
  double m_conv_force;
  double m_conv_velocity;

public:
  LBFluidNode() {
    add_parameters(
        {{"_index", AutoParameter::read_only, [this]() { return m_index; }}});
  }

  void do_construct(VariantMap const &params) override {
    auto const lb_sip =
        get_value<std::shared_ptr<LBFluid>>(params, "parent_sip");
    m_lb_fluid = lb_sip->get_lb_fluid();
    auto const lb_params = lb_sip->get_lb_params();
    auto const tau = lb_params->get_tau();
    auto const agrid = lb_params->get_agrid();
    m_conv_dens = Utils::int_pow<3>(agrid);
    m_conv_press = Utils::int_pow<1>(agrid) * Utils::int_pow<2>(tau);
    m_conv_force = Utils::int_pow<2>(tau) / Utils::int_pow<1>(agrid);
    m_conv_velocity = Utils::int_pow<1>(tau) / Utils::int_pow<1>(agrid);
    m_grid_size = m_lb_fluid->get_lattice().get_grid_dimensions();
    m_index = get_mapped_index(get_value<Utils::Vector3i>(params, "index"),
                               m_grid_size);
  }

  Variant do_call_method(std::string const &name,
                         VariantMap const &params) override;
};
} // namespace ScriptInterface::walberla

#endif // WALBERLA
