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
#ifndef UTILS_TYPE_TRAITS_HPP
#define UTILS_TYPE_TRAITS_HPP

#include <climits>
#include <type_traits>

namespace Utils {

/**
 * @brief Remove const from a function signature.
 */
template <typename T> struct function_remove_const;

template <typename R, typename... Args>
struct function_remove_const<R(Args...)> {
  using type = R(Args...);
};

template <typename R, typename... Args>
struct function_remove_const<R(Args...) const> {
  using type = R(Args...);
};

template <class...> struct conjunction : std::true_type {};
template <class B1> struct conjunction<B1> : B1 {};
template <class B1, class... Bn>
struct conjunction<B1, Bn...>
    : std::conditional<bool(B1::value), conjunction<Bn...>, B1>::type {};

/**
 * @brief Size of a type in bits.
 *
 * On posix platforms this is 8 * sizeof(T).
 */
template <class T>
struct size_in_bits : std::integral_constant<size_t, CHAR_BIT * sizeof(T)> {};
} // namespace Utils

#endif
