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
#include <vector>

namespace fft {
namespace detail {
void fft_free(void *p);
void *fft_malloc(std::size_t length);
} // namespace detail

/** @brief Aligned allocator for FFT data. */
template <class T> struct allocator {
  typedef T value_type;
  allocator() noexcept = default; // default ctor not required
  template <class U> explicit allocator(const allocator<U> &) {}
  template <class U> bool operator==(const allocator<U> &) const {
    return true;
  }
  template <class U> bool operator!=(const allocator<U> &) const {
    return false;
  }

  T *allocate(const std::size_t n) const {
    if (n == 0) {
      return nullptr;
    }
    if (n > std::numeric_limits<std::size_t>::max() / sizeof(T)) {
      throw std::bad_array_new_length();
    }
    void *const pv = detail::fft_malloc(n * sizeof(T));
    if (!pv) {
      throw std::bad_alloc();
    }
    return static_cast<T *>(pv);
  }

  void deallocate(T *const p, std::size_t) const noexcept {
    detail::fft_free(static_cast<void *>(p));
  }
};

template <class T> using vector = std::vector<T, allocator<T>>;

} // namespace fft
