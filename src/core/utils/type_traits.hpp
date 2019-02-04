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

template <bool P, typename T = void>
using enable_if_t = typename std::enable_if<P, T>::type;

/**
 * @brief Implementation of std::void_t
 *
 * See @url https://en.cppreference.com/w/cpp/types/void_t
 *
 */
template <class... Ts>
using void_t = void;

    namespace detail {
        template <template <class...> class Trait, class Enabler, class... Args>
        struct is_detected : std::false_type{};

        template <template <class...> class Trait, class... Args>
        struct is_detected<Trait, void_t<Trait<Args...>>, Args...> : std::true_type{};
    }

    /**
     * @brief Partial implementation of std::experimental::is_detected.
     *
     * Source and in-depth discussion at
     * @url https://blog.tartanllama.xyz/detection-idiom and
     * @url https://en.cppreference.com/w/cpp/experimental/is_detected
     *
     * This can be used to check if an expression is well-formed for a
     * set of types, like so:
     *
     * template<typename T>
     * using get_t = decltype(std::declval<T>().get());
     * template<typename T>
     * using has_get = is_detected<get_t, T>;
     *
     */
    template <template <class...> class Trait, class... Args>
    using is_detected = typename detail::is_detected<Trait, void, Args...>::type;

} // namespace Utils

#endif
