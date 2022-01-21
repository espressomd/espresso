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
#ifndef SCRIPT_INTERFACE_WALBERLA_OPTIONAL_REDUCTION_HPP
#define SCRIPT_INTERFACE_WALBERLA_OPTIONAL_REDUCTION_HPP

#include "config.hpp"

#ifdef LB_WALBERLA

#include "script_interface/None.hpp"
#include "script_interface/Variant.hpp"

#include "core/communication.hpp"

#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/optional.hpp>
#include <boost/serialization/vector.hpp>

#include <cassert>
#include <functional>
#include <type_traits>
#include <utility>

namespace ScriptInterface::walberla {

/**
 * @brief Reduction of an optional value with unit conversion.
 *
 * @param result       optional value to reduce on main rank
 * @param conversion   unit conversion
 *
 * @return Reduced optional wrapped in a variant (on head node) or empty
 *         variant (on worker nodes)
 */
template <typename T, bool Conversion = true>
Variant optional_reduction_with_conversion(boost::optional<T> const &result,
                                           double conversion) {
  assert(1 == boost::mpi::all_reduce(comm_cart, static_cast<int>(!!result),
                                     std::plus<>()) &&
         "Incorrect number of return values");
  auto const conversion_handler = [=](T const &value) {
    if constexpr (Conversion) {
      return Variant{value * (1. / conversion)};
    } else {
      return Variant{value};
    }
  };
  if (comm_cart.rank() == 0) {
    if (result) {
      return conversion_handler(*result);
    }
    T value;
    comm_cart.recv(boost::mpi::any_source, 42, value);
    return conversion_handler(value);
  } else if (result) {
    comm_cart.send(0, 42, *result);
  }
  return {};
}

template <typename T>
Variant optional_reduction_with_conversion(boost::optional<T> const &result) {
  return optional_reduction_with_conversion<T, false>(result, -1.);
}

} // namespace ScriptInterface::walberla

#endif // LB_WALBERLA
#endif
