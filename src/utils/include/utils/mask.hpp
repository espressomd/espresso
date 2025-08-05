/*
 * Copyright (C) 2019-2022 The ESPResSo project
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

#include <concepts>
#include <cstddef>
#include <limits>
#include <type_traits>
#include <utility>

namespace Utils {
namespace detail {
template <class T, std::size_t... I>
auto mask_impl(std::unsigned_integral auto mask, T const &t,
               std::index_sequence<I...>) {
  return T{((mask & (1u << I)) ? get<I>(t) : std::tuple_element_t<I, T>{})...};
}
} // namespace detail

/**
 * @brief Pick elements of a tuple-like by a bit mask.
 *
 * E.g. every element of the input for which the corresponding
 * bit is set in the mask is set is copied to the output unmodified,
 * the elements that are not set are set to zero (default constructed
 * instance of the type).
 *
 * Example:
 *   <tt>mask(0b1011u, {1, 2, 3, 4}) => {1, 2, 0, 4}</tt>
 *
 * @tparam T implements the tuple interface(get, tuple_size, ...)
 * @param mask bit mask, if the i-th bit is set, the i-th element
 *        in @p t is copied to the output, otherwise it is set to zero.
 * @param t input elements
 * @return t partially zeroed out according to mask
 */
template <class T> T mask(std::unsigned_integral auto mask, T const &t) {
  auto constexpr size_in_bits = std::numeric_limits<decltype(mask)>::digits;
  static_assert(size_in_bits >= std::tuple_size_v<T>);
  return detail::mask_impl(mask, t,
                           std::make_index_sequence<std::tuple_size_v<T>>{});
}
} // namespace Utils
