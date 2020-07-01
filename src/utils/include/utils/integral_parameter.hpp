/*
 * Copyright (C) 2020 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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
#ifndef ESPRESSO_INTEGRAL_PARAMETER_HPP
#define ESPRESSO_INTEGRAL_PARAMETER_HPP

#include <stdexcept>

namespace Utils {
namespace detail {
template <template <size_t> class F, size_t I, size_t N>
struct integral_parameter_impl {
  template <class... Args>
  static decltype(auto) eval(size_t i, Args &&... args) {
    if (i == I)
      return F<I>{}(std::forward<Args>(args)...);

    return integral_parameter_impl<F, I + 1, N>::eval(
        i, std::forward<Args>(args)...);
  }
};

template <template <size_t> class F, size_t N>
struct integral_parameter_impl<F, N, N> {
  template <class... Args>
  static decltype(auto) eval(size_t i, Args &&... args) {
    if (i == N)
      return F<N>{}(std::forward<Args>(args)...);

    throw std::runtime_error("Invalid parameter value");
  }
};
} // namespace detail

/**
 * @brief Generate a call table for a integral non-type template parameter.
 */
template <template <size_t> class F, size_t min, size_t max, class... Args>
decltype(auto) integral_parameter(size_t i, Args &&... args) {
  return detail::integral_parameter_impl<F, min, max>::eval(
      i, std::forward<Args>(args)...);
}

} // namespace Utils

#endif // ESPRESSO_INTEGRAL_PARAMETER_HPP
