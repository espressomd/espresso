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

#include "LBFluidSlice.hpp"

#include "LatticeSlice.impl.hpp"

#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace ScriptInterface::walberla {

Variant LBFluidSlice::do_call_method(std::string const &name,
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
  if (name == "get_lattice_speed") {
    return 1. / m_conv_velocity;
  }

  // slice getter/setter callback
  auto const call = [this, &params](auto method_ptr,
                                    std::vector<int> const &data_dims,
                                    double units = 1.) -> Variant {
    auto &obj = *m_lb_fluid;
    if constexpr (std::is_invocable_v<decltype(method_ptr), LatticeModel *,
                                      Utils::Vector3i const &,
                                      Utils::Vector3i const &>) {
      return gather_3d(params, data_dims, obj, method_ptr, units);
    } else {
      scatter_3d(params, data_dims, obj, method_ptr, units);
      return {};
    }
  };

  if (name == "get_population") {
    auto const pop_size = m_shape_val.at("population");
    return call(&LatticeModel::get_slice_population, pop_size);
  }
  if (name == "set_population") {
    auto const pop_size = m_shape_val.at("population");
    return call(&LatticeModel::set_slice_population, pop_size);
  }
  if (name == "get_density") {
    return call(&LatticeModel::get_slice_density, {1}, 1. / m_conv_dens);
  }
  if (name == "set_density") {
    return call(&LatticeModel::set_slice_density, {1}, m_conv_dens);
  }
  if (name == "get_velocity") {
    return call(&LatticeModel::get_slice_velocity, {3}, 1. / m_conv_velocity);
  }
  if (name == "set_velocity") {
    return call(&LatticeModel::set_slice_velocity, {3}, m_conv_velocity);
  }
  if (name == "get_is_boundary") {
    return call(&LatticeModel::get_slice_is_boundary, {1});
  }
  if (name == "get_velocity_at_boundary") {
    return call(&LatticeModel::get_slice_velocity_at_boundary, {1},
                1. / m_conv_velocity);
  }
  if (name == "set_velocity_at_boundary") {
    return call(&LatticeModel::set_slice_velocity_at_boundary, {1},
                m_conv_velocity);
  }
  if (name == "get_pressure_tensor") {
    return call(&LatticeModel::get_slice_pressure_tensor, {3, 3},
                1. / m_conv_press);
  }
  if (name == "get_pressure_tensor_neq") {
    auto variant = do_call_method("get_pressure_tensor", params);
    if (context()->is_head_node()) {
      auto constexpr c_sound_sq = 1. / 3.;
      auto const density = m_lb_fluid->get_density();
      auto const diagonal_term = density * c_sound_sq / m_conv_press;
      // modify existing variant in-place
      auto &vec = *(boost::get<std::vector<double>>(
          &(boost::get<std::vector<Variant>>(&variant)->at(0))));
      for (auto it = vec.begin(); it < vec.end(); it += 9) {
        *(it + 0) -= diagonal_term;
        *(it + 4) -= diagonal_term;
        *(it + 8) -= diagonal_term;
      }
    }
    return variant;
  }
  if (name == "get_last_applied_force") {
    return call(&LatticeModel::get_slice_last_applied_force, {3},
                1. / m_conv_force);
  }
  if (name == "set_last_applied_force") {
    return call(&LatticeModel::set_slice_last_applied_force, {3}, m_conv_force);
  }

  return {};
}

} // namespace ScriptInterface::walberla

#endif // WALBERLA
