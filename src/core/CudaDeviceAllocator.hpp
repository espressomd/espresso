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
#ifndef ESPRESSO_CUDA_DEVICE_ALLOCATOR_HPP
#define ESPRESSO_CUDA_DEVICE_ALLOCATOR_HPP

#include <utils/device_qualifier.hpp>

#ifndef __CUDACC__
#error "This header can only be used with a CUDA-capable compiler."
#endif

#include <thrust/device_free.h>
#include <thrust/device_malloc.h>
#include <thrust/device_ptr.h>
#include <thrust/device_reference.h>

#include <cstddef>

/**
 * @brief Allocator that uses CUDA to allocate GPU memory.
 *
 * This allocator uses thrust::malloc and thrust::free to
 * manage device memory. To be able to have static containers
 * for device memory, exceptions in deallocations are suppressed.
 * This is needed because static containers may be destroyed after
 * the CUDA system was teared down, leading to device_free to fail.
 *
 * @tparam T Type to allocate memory for.
 */
template <class T> struct CudaDeviceAllocator {
  using value_type = T;
  using pointer = thrust::device_ptr<T>;
  using const_pointer = thrust::device_ptr<const T>;
  using reference = thrust::device_reference<T>;
  using const_reference = thrust::device_reference<const T>;

  CudaDeviceAllocator() noexcept = default;
  template <class U>
  explicit CudaDeviceAllocator(const CudaDeviceAllocator<U> &) {}
  template <class U> bool operator==(const CudaDeviceAllocator<U> &) const {
    return true;
  }
  template <class U> bool operator!=(const CudaDeviceAllocator<U> &) const {
    return false;
  }

  pointer allocate(const std::size_t n) const {
    return thrust::device_malloc<T>(n);
  }
  void deallocate(pointer p, std::size_t) const noexcept {
    try {
      thrust::device_free(p);
    } catch (thrust::system::system_error const &) {
      ;
    }
  }
};
#endif
