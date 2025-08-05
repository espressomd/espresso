/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#include <iterator>
#include <type_traits>

namespace Utils {
namespace detail {
template <class ContainerOrValue, class OutputIterator>
OutputIterator flatten_impl(ContainerOrValue const &c, OutputIterator out) {
  if constexpr (std::is_assignable_v<decltype(*out), ContainerOrValue>) {
    *out = c;
    return ++out;
  } else {
    using ValueType = typename ContainerOrValue::value_type;
    for (auto const &e : c) {
      out = flatten_impl<ValueType, OutputIterator>(e, out);
    }
    return out;
  }
}
} // namespace detail

/**
 * @brief Flatten a range of ranges.
 *
 * Copy a range of ranges to an output range by subsequently
 * copying the nested ranges to the output. Arbitrary deep
 * nesting is supported, the elements are copied into the output
 * in a depth-first fashion.
 *
 * @tparam Range A Forward Range
 * @tparam OutputIterator An OutputIterator
 * @param v Input Range
 * @param out Output iterator
 */
template <class Range, class OutputIterator>
void flatten(Range const &v, OutputIterator out) {
  detail::flatten_impl<Range, OutputIterator>(v, out);
}
} // namespace Utils
