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

#include "EKSpeciesSlice.hpp"

#include "LatticeSlice.impl.hpp"

#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace ScriptInterface::walberla {

Variant EKSpeciesSlice::do_call_method(std::string const &name,
                                       VariantMap const &params) {
  if (name == "get_slice_size") {
    return {m_slice_upper_corner - m_slice_lower_corner};
  }
  if (name == "get_slice_ranges") {
    return {std::vector<Variant>{m_slice_lower_corner, m_slice_upper_corner}};
  }
  if (name == "get_ek_sip") {
    return {m_ek_sip};
  }
  if (name == "get_value_shape") {
    auto const name = get_value<std::string>(params, "name");
    if (m_shape_val.count(name) == 0) {
      context()->parallel_try_catch([&]() {
        throw std::runtime_error("Unknown EK property '" + name + "'");
      });
    }
    return m_shape_val.at(name);
  }

  // slice getter/setter callback
  auto const call = [this, params](auto method_ptr,
                                   std::vector<int> const &data_dims,
                                   double units = 1.) -> Variant {
    auto &obj = *m_ek_species;
    if constexpr (std::is_invocable_v<decltype(method_ptr), LatticeModel *,
                                      Utils::Vector3i const &,
                                      Utils::Vector3i const &>) {
      return gather_3d(params, data_dims, obj, method_ptr, units);
    } else {
      scatter_3d(params, data_dims, obj, method_ptr, units);
      return {};
    }
  };

  if (name == "get_density") {
    return call(&LatticeModel::get_slice_density, {1}, 1. / m_conv_dens);
  }
  if (name == "set_density") {
    return call(&LatticeModel::set_slice_density, {1}, m_conv_dens);
  }
  if (name == "get_is_boundary") {
    return call(&LatticeModel::get_slice_is_boundary, {1});
  }
  if (name == "get_flux_at_boundary") {
    return call(&LatticeModel::get_slice_flux_at_boundary, {1},
                1. / m_conv_flux);
  }
  if (name == "set_flux_at_boundary") {
    return call(&LatticeModel::set_slice_flux_boundary, {1}, m_conv_flux);
  }
  if (name == "get_density_at_boundary") {
    return call(&LatticeModel::get_slice_density_at_boundary, {1},
                1. / m_conv_dens);
  }
  if (name == "set_density_at_boundary") {
    return call(&LatticeModel::set_slice_density_boundary, {1}, m_conv_dens);
  }

  return {};
}

} // namespace ScriptInterface::walberla

#endif // WALBERLA
