/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016,2017,2018 The ESPResSo
  project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
  Max-Planck-Institute for Polymer Research, Theory Group

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

#ifndef UTILS_MATH_INT_POW_HPP
#define UTILS_MATH_INT_POW_HPP

#include "utils/device_qualifier.hpp"

#include <type_traits>

namespace Utils {
namespace detail {
template <class T, unsigned n, class = void> struct int_pow_impl {
  DEVICE_QUALIFIER constexpr T operator()(T x) const {
    return x * int_pow_impl<T, (n - 1) / 2>{}(x * x);
  }
};

/* Specializaton for n even */
template <class T, unsigned n>
struct int_pow_impl<T, n, std::enable_if_t<n % 2 == 0>> {
  DEVICE_QUALIFIER constexpr T operator()(T x) const {
    return int_pow_impl<T, n / 2>{}(x * x);
  }
};

template <class T> struct int_pow_impl<T, 1> {
  DEVICE_QUALIFIER constexpr T operator()(T x) const { return x; }
};

template <class T> struct int_pow_impl<T, 0> {
  DEVICE_QUALIFIER constexpr T operator()(T x) const { return T{1}; }
};
} // namespace detail

/**
 * \brief Calculate integer powers.
 * This functions calculates x^n, where
 * n is a positive integer that is known
 * at compile time. It uses exponentiation by
 * squaring to construct a efficient function.
 */
template <unsigned n, typename T> DEVICE_QUALIFIER constexpr T int_pow(T x) {
  return detail::int_pow_impl<T, n>{}(x);
}
} // namespace Utils

#endif
