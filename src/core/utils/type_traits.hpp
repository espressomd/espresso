/*
Copyright (C) 2010-2018 The ESPResSo project

This file is part of ESPResSo.

ESPResSo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ESPResSo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef UTILS_TYPE_TRAITS_HPP
#define UTILS_TYPE_TRAITS_HPP

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

template <class T>
using function_remove_const_t = typename function_remove_const<T>::type;

/**
 * @brief True iff T is an instantiation of Template.
 */
template <typename T, template <typename...> class Template>
struct is_instance_of : public std::false_type {};

template <typename... T, template <typename...> class Template>
struct is_instance_of<Template<T...>, Template> : public std::true_type {};

template <bool P, typename T = void>
using enable_if_t = typename std::enable_if<P, T>::type;

template <typename T> using decay_t = typename std::decay<T>::type;

template <class T> using add_const_t = typename std::add_const<T>::type;

template <class...> struct conjunction : std::true_type {};
template <class B1> struct conjunction<B1> : B1 {};
template <class B1, class... Bn>
struct conjunction<B1, Bn...>
    : std::conditional<bool(B1::value), conjunction<Bn...>, B1>::type {};

} // namespace Utils

#endif
