/*
 * Copyright (C) 2010-2022 The ESPResSo project
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
#include <type_traits>
#include <vector>

void cuda_malloc_host(void **p, std::size_t n);
void cuda_free_host(void *p);

/**
 * @brief Allocator that uses CUDA to allocate CPU memory.
 *
 * Using the CUDA allocator can have performance benefits,
 * because it returns pinned memory that is suitable for
 * DMA.
 *
 * @tparam T Type to allocate memory for.
 */
template <class T> struct CudaHostAllocator {
  using value_type = T;
  using pointer = T *;
  using reference = T &;
  using const_reference = std::add_const_t<reference>;

  CudaHostAllocator() noexcept = default;
  template <class U> explicit CudaHostAllocator(const CudaHostAllocator<U> &) {}
  template <class U> bool operator==(const CudaHostAllocator<U> &) const {
    return true;
  }
  template <class U> bool operator!=(const CudaHostAllocator<U> &) const {
    return false;
  }

  T *allocate(const std::size_t n) const {
    T *result{};

    cuda_malloc_host(reinterpret_cast<void **>(&result),
                     n * sizeof(value_type));

    return result;
  }
  void deallocate(T *const p, std::size_t) const noexcept { cuda_free_host(p); }
};

template <class T> using pinned_vector = std::vector<T, CudaHostAllocator<T>>;
