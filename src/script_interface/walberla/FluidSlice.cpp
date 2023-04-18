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

#include "config/config.hpp"

#ifdef WALBERLA

#include "FluidSlice.hpp"

#include "LatticeSlice.impl.hpp"

#include <walberla_bridge/lattice_boltzmann/LBWalberlaBase.hpp>

#include <stdexcept>
#include <string>
#include <vector>

namespace ScriptInterface::walberla {

Variant FluidSlice::do_call_method(std::string const &name,
                                   VariantMap const &params) {
  if (name == "get_slice_size") {
    return {m_slice_upper_corner - m_slice_lower_corner};
  }
  if (name == "get_slice_ranges") {
    return {std::vector<Variant>{m_slice_lower_corner, m_slice_upper_corner}};
  }
  if (name == "get_lb_sip") {
    return {m_lb_sip};
  }
  if (name == "get_value_shape") {
    auto const name = get_value<std::string>(params, "name");
    if (m_shape_val.count(name) == 0) {
      context()->parallel_try_catch([&]() {
        throw std::runtime_error("Unknown fluid property '" + name + "'");
      });
    }
    return m_shape_val.at(name);
  }
  if (name == "get_density") {
    return gather_3d(params, {1}, *m_lb_fluid,
                     &::LBWalberlaBase::get_slice_density, 1. / m_conv_dens);
  }
  if (name == "set_density") {
    scatter_3d(params, {1}, *m_lb_fluid, &::LBWalberlaBase::set_slice_density,
               m_conv_dens);
    return {};
  }
  if (name == "get_population") {
    return gather_3d(params, m_shape_val.at("population"), *m_lb_fluid,
                     &::LBWalberlaBase::get_slice_population);
  }
  if (name == "set_population") {
    scatter_3d(params, m_shape_val.at("population"), *m_lb_fluid,
               &::LBWalberlaBase::set_slice_population);
    return {};
  }
  if (name == "get_velocity") {
    return gather_3d(params, {3}, *m_lb_fluid,
                     &::LBWalberlaBase::get_slice_velocity,
                     1. / m_conv_velocity);
  }
  if (name == "set_velocity") {
    scatter_3d(params, {3}, *m_lb_fluid, &::LBWalberlaBase::set_slice_velocity,
               m_conv_velocity);
    return {};
  }
  if (name == "get_is_boundary") {
    return gather_3d(params, {1}, *m_lb_fluid,
                     &::LBWalberlaBase::get_slice_is_boundary);
  }
  if (name == "get_velocity_at_boundary") {
    return gather_3d(params, {1}, *m_lb_fluid,
                     &::LBWalberlaBase::get_slice_velocity_at_boundary,
                     1. / m_conv_velocity);
  }
  if (name == "set_velocity_at_boundary") {
    scatter_3d(params, {1}, *m_lb_fluid,
               &::LBWalberlaBase::set_slice_velocity_at_boundary,
               m_conv_velocity);
    return {};
  }
  if (name == "get_pressure_tensor") {
    return gather_3d(params, {3, 3}, *m_lb_fluid,
                     &::LBWalberlaBase::get_slice_pressure_tensor,
                     1. / m_conv_press);
  }
  if (name == "get_last_applied_force") {
    return gather_3d(params, {3}, *m_lb_fluid,
                     &::LBWalberlaBase::get_slice_last_applied_force,
                     1. / m_conv_force);
  }
  if (name == "set_last_applied_force") {
    scatter_3d(params, {3}, *m_lb_fluid,
               &::LBWalberlaBase::set_slice_last_applied_force, m_conv_force);
    return {};
  }
  if (name == "get_lattice_speed") {
    return 1. / m_conv_velocity;
  }

  return {};
}

} // namespace ScriptInterface::walberla

#endif // WALBERLA
