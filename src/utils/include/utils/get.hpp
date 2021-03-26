/*
 * Copyright (C) 2010-2019 The ESPResSo project
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
#ifndef UTILS_GET_HPP
#define UTILS_GET_HPP

#include <cstddef>
#include <tuple>

namespace Utils {
template <std::size_t I, typename T>
const std::tuple_element_t<I, T> &get(const T &v) {
  return std::get<I>(v);
}

template <class T> struct tuple_size : std::tuple_size<T> {};

template <std::size_t I, class Tuple>
struct tuple_element : std::tuple_element<I, Tuple> {};

template <std::size_t I, class Tuple>
using tuple_element_t = typename tuple_element<I, Tuple>::type;
} // namespace Utils
#endif
