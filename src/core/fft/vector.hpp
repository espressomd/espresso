/*
 * Copyright (C) 2010-2024 The ESPResSo project
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

#pragma once

#include <cstddef>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <vector>

namespace fft {

/** @brief Aligned allocator for FFT data. */
template <class FloatType> struct allocator {
  static_assert(std::is_same_v<FloatType, float> or
                    std::is_same_v<FloatType, double>,
                "FFTW only implements float and double");
  typedef FloatType value_type;
  allocator() noexcept = default; // default ctor not required
  template <class U> explicit allocator(allocator<U> const &) {}
  template <class U> bool operator==(allocator<U> const &) const {
    return true;
  }
  template <class U> bool operator!=(allocator<U> const &) const {
    return false;
  }

  FloatType *allocate(std::size_t n) const;

  void deallocate(FloatType *p, std::size_t) const noexcept;
};

template <class T> using vector = std::vector<T, allocator<T>>;

} // namespace fft
