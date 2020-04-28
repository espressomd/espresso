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
      (void)f(e);
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
 * @tparam Pred Must be callable with all types occurring as elements
 *              in Tuple, returning a value convertible to bool.
 * @tparam Tuple a tuple-like type supporting Utils::get
 * @tparam F Must be callable with all types occurring as elements
 *           in Tuple.
 * @param pred The predicate.
 * @param t Tuple for search in.
 * @param f Is called for first found element.
 * @return true iff the predicate was true for at least one element.
 */
template <class Pred, class Tuple, class F>
constexpr auto find_if(Pred &&pred, Tuple const &t, F &&f) {
  return detail::find_if_impl<0, Utils::tuple_size<Tuple>::value>::eval(
      std::forward<Pred>(pred), t, std::forward<F>(f));
}

namespace detail {
template <template <class> class Predicate, size_t I, size_t N>
struct filter_impl {
  template <class Tuple>
  constexpr static auto get(Tuple const &t, std::true_type) {
    using Utils::get;

    return std::make_tuple(get<I>(t));
  }

  template <class Tuple>
  constexpr static auto get(Tuple const &t, std::false_type) {
    return std::make_tuple();
  }

  template <class Tuple> constexpr static auto eval(Tuple const &t) {
    using Utils::tuple_element_t;
    using element_type = tuple_element_t<I, Tuple>;
    constexpr bool pred = Predicate<element_type>::value;

    return std::tuple_cat(
        get(t, std::conditional_t<pred, std::true_type, std::false_type>{}),
        filter_impl<Predicate, I + 1, N>::eval(t));
  }
};

template <template <class> class Predicate, size_t I>
struct filter_impl<Predicate, I, I> {
  template <class Tuple> constexpr static auto eval(Tuple const &) {
    return std::make_tuple();
  }
};

} // namespace detail

/**
 * @brief Filter a tuple by a static predicate.
 *
 * E.g. filter(std::is_integral, std::make_tuple(1, 1.5, 2u))
 *       -> std::tuple<int, unsigned>(1, 2u).
 *
 *       @tparam Predicate Metafunction that returns a bool for each type.
 *       @tparam Tuple tuple-like type that implements the tuple interface.
 *       @param in The input tuple-like
 *
 *       @return std::tuple with the elements selected by the predicate.
 */
template <template <class> class Predicate, class Tuple>
constexpr auto filter(Tuple const &in) {
  using Utils::tuple_size;

  return detail::filter_impl<Predicate, 0, tuple_size<Tuple>::value>::eval(in);
}
} // namespace Utils

#endif // ESPRESSO_TUPLE_HPP
