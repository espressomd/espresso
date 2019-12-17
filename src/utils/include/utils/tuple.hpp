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

#include <utils/get.hpp>
#include <utils/type_traits.hpp>

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

namespace detail {
template <size_t I, size_t N> struct find_if_impl {
  template <class Pred, class Tuple, class F>
  static constexpr bool eval(Pred &&pred, Tuple const &t, F &&f) {
    using Utils::get;
    auto const &e = get<I>(t);

    if (pred(e)) {
      f(e);
      return true;
    }

    return find_if_impl<I + 1, N>::eval(std::forward<Pred>(pred), t,
                                        std::forward<F>(f));
  }
};

template <size_t N> struct find_if_impl<N, N> {
  template <class Pred, class Tuple, class F>
  static constexpr bool eval(Pred, Tuple, F) {
    return false;
  }
};
} // namespace detail

/**
 * @brief Find first element in tuple for which a predicate holds,
 *        and call a Callable with the element.
 *
 * Any return value of the Callable is ignored.
 *
 * @tparam Pred Must be callable with all type occurring as elements
 *              in Tuple, returning a value convertible to bool.
 * @tparam Tuple a tuple-like type supporting Utils::get
 * @tparam F Must be callable with all type occurring as elements
 *           in Tuple.
 * @param pred The predicate.
 * @param t Tuple for search in.
 * @param f Is called for first found element.
 * @return Result of invocation of f if an element was found,
 *         nothing otherwise.
 */
template <class Pred, class Tuple, class F>
constexpr auto find_if(Pred &&pred, Tuple const &t, F &&f) {
  return detail::find_if_impl<0, Utils::tuple_size<Tuple>::value>::eval(
      std::forward<Pred>(pred), t, std::forward<F>(f));
}
} // namespace Utils

#endif // ESPRESSO_TUPLE_HPP
