/*
 * Copyright (C) 2024-2026 The ESPResSo project
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

#include <config/config.hpp>

#ifdef ESPRESSO_WALBERLA

#include "EKIndexedReactionSlice.hpp"

#include "LatticeSlice.impl.hpp"

#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace ScriptInterface::walberla {

Variant EKIndexedReactionSlice::do_call_method(std::string const &name,
                                               VariantMap const &params) {
  if (name == "get_slice_size") {
    return {m_slice_upper_corner - m_slice_lower_corner};
  }
  if (name == "get_slice_ranges") {
    return {std::vector<Variant>{m_slice_lower_corner, m_slice_upper_corner}};
  }
  if (name == "get_reaction_sip") {
    return {m_reaction_sip};
  }
  if (name == "get_value_shape") {
    auto const attr = get_value<std::string>(params, "name");
    if (not m_shape_val.contains(attr)) {
      context()->parallel_try_catch([&]() {
        throw std::runtime_error("Unknown EK indexed reaction property '" +
                                 attr + "'");
      });
    }
    return m_shape_val.at(attr);
  }

  auto adapter = EKReactionSliceAdapter{m_reaction_impl.get()};

  auto const call = [this, &adapter,
                     &params](auto method_ptr,
                              std::vector<int> const &data_dims) -> Variant {
    if constexpr (std::is_invocable_v<
                      decltype(method_ptr), EKReactionSliceAdapter *,
                      Utils::Vector3i const &, Utils::Vector3i const &>) {
      return gather_3d(data_dims, adapter, method_ptr);
    } else {
      scatter_3d(params.at("values"), data_dims, adapter, method_ptr);
      return {};
    }
  };

  if (name == "get_is_boundary") {
    return call(&EKReactionSliceAdapter::get_slice_is_boundary, {1});
  }
  if (name == "set_is_boundary") {
    return call(&EKReactionSliceAdapter::set_slice_is_boundary, {1});
  }

  return {};
}

} // namespace ScriptInterface::walberla

#endif // ESPRESSO_WALBERLA
