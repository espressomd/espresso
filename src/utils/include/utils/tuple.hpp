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
#ifndef ESPRESSO_TUPLE_HPP
#define ESPRESSO_TUPLE_HPP

/**
 * @file
 * Algorithms for tuple-like inhomogeneous containers.
 */

namespace Utils {
namespace detail {
template <class F, class Tuple, std::size_t... I>
constexpr decltype(auto) apply_impl(F &&f, Tuple &&t,
                                    std::index_sequence<I...>) {
  return f(std::get<I>(std::forward<Tuple>(t))...);
}
} // namespace detail

/**
 * @brief Call function with expanded tuple as parameters.
 *
 * Like std::apply.
 *
 * @tparam F Callable with tuple elements as arguments
 * @tparam Tuple Has to conform to the tuple interface
 * @return Whatever @p f returns.
 */
template <class F, class Tuple>
constexpr decltype(auto) apply(F &&f, Tuple &&t) {
  return detail::apply_impl(
      std::forward<F>(f), std::forward<Tuple>(t),
      std::make_index_sequence<
          std::tuple_size<std::remove_reference_t<Tuple>>::value>{});
}

namespace detail {
template <class Tuple, class F, std::size_t... I>
constexpr void for_each_impl(F &&f, Tuple t, std::index_sequence<I...>) {
  using expand = int[];
  (void)expand{0, ((void)(f(std::get<I>(t))), 0)...};
}

template <class F, class Tuple>
constexpr void for_each_impl(F, Tuple, std::index_sequence<>) {}
} // namespace detail

template <class F, class Tuple> void for_each(F &&f, Tuple &t) {
  detail::for_each_impl<Tuple &>(
      std::forward<F>(f), t,
      std::make_index_sequence<std::tuple_size<Tuple>::value>{});
}

template <class F, class Tuple> void for_each(F &&f, Tuple &&t) {
  detail::for_each_impl(
      std::forward<F>(f), std::forward<Tuple>(t),
      std::make_index_sequence<
          std::tuple_size<std::remove_reference_t<Tuple>>::value>{});
}
} // namespace Utils

#endif // ESPRESSO_TUPLE_HPP
