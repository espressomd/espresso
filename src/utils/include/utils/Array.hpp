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
#ifndef UTILS_ARRAY_HPP
#define UTILS_ARRAY_HPP

#include "device_qualifier.hpp"
#include "get.hpp"

#include <boost/serialization/access.hpp>

#include <cassert>
#include <cstddef>
#include <stdexcept>

namespace Utils {

namespace detail {

template <typename T, std::size_t N> struct Storage {
  T m_data[N];

private:
  friend boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, const unsigned int /* version */) {
    ar &m_data;
  }
};

template <typename T> struct Storage<T, 0> {};

} // namespace detail

template <typename T, std::size_t N> struct Array {
  using value_type = T;
  using size_type = std::size_t;
  using difference_type = std::ptrdiff_t;
  using reference = value_type &;
  using const_reference = const value_type &;
  using iterator = value_type *;
  using const_iterator = const value_type *;
  using pointer = value_type *;
  using const_pointer = const value_type *;

  detail::Storage<T, N> m_storage;

  DEVICE_QUALIFIER constexpr reference at(size_type i) {
    if (i >= N) {
      DEVICE_THROW(std::out_of_range("Array access out of bounds."));
    }
    return m_storage.m_data[i];
  }

  DEVICE_QUALIFIER constexpr const_reference at(size_type i) const {
    if (i >= N) {
      DEVICE_THROW(std::out_of_range("Array access out of bounds."));
    }
    return m_storage.m_data[i];
  }

  DEVICE_QUALIFIER constexpr reference operator[](size_type i) {
    DEVICE_ASSERT(i < N);
    return m_storage.m_data[i];
  }

  DEVICE_QUALIFIER constexpr const_reference operator[](size_type i) const {
    DEVICE_ASSERT(i < N);
    return m_storage.m_data[i];
  }

  DEVICE_QUALIFIER constexpr reference front() { return *begin(); }

  DEVICE_QUALIFIER constexpr const_reference front() const { return *cbegin(); }

  DEVICE_QUALIFIER constexpr reference back() { return *(end() - 1); }

  DEVICE_QUALIFIER constexpr const_reference back() const {
    return *(cend() - 1);
  }

  DEVICE_QUALIFIER constexpr pointer data() noexcept {
    return &m_storage.m_data[0];
  }

  DEVICE_QUALIFIER constexpr const_pointer data() const noexcept {
    return &m_storage.m_data[0];
  }

  DEVICE_QUALIFIER constexpr iterator begin() noexcept {
    return &m_storage.m_data[0];
  };

  DEVICE_QUALIFIER constexpr const_iterator begin() const noexcept {
    return &m_storage.m_data[0];
  };

  DEVICE_QUALIFIER constexpr const_iterator cbegin() const noexcept {
    return &m_storage.m_data[0];
  };

  DEVICE_QUALIFIER constexpr iterator end() noexcept {
    return &m_storage.m_data[N];
  };

  DEVICE_QUALIFIER constexpr const_iterator end() const noexcept {
    return &m_storage.m_data[N];
  };

  DEVICE_QUALIFIER constexpr const_iterator cend() const noexcept {
    return &m_storage.m_data[N];
  };

  DEVICE_QUALIFIER constexpr bool empty() const noexcept { return size() == 0; }

  DEVICE_QUALIFIER constexpr size_type size() const noexcept { return N; }

  DEVICE_QUALIFIER constexpr size_type max_size() const noexcept { return N; }

  DEVICE_QUALIFIER void fill(const value_type &value) {
    for (size_type i = 0; i < size(); ++i)
      m_storage.m_data[i] = value;
  }

  DEVICE_QUALIFIER static constexpr Array<T, N>
  broadcast(const value_type &value) {
    Array<T, N> ret{};
    for (size_type i = 0; i < N; ++i) {
      ret[i] = value;
    }
    return ret;
  }

#ifdef __HCC__
  // workaround for https://github.com/RadeonOpenCompute/hcc/issues/860
  __attribute__((annotate("serialize"))) void
  __cxxamp_serialize(Kalmar::Serialize &s) const;
  __attribute__((annotate("user_deserialize"))) void cxxamp_deserialize()
      [[cpu]] [[hc]];
#endif

private:
  friend boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, const unsigned int /* version */) {
    ar &m_storage;
  }
};

template <std::size_t I, class T, std::size_t N>
struct tuple_element<I, Array<T, N>> {
  using type = T;
};

template <class T, std::size_t N>
struct tuple_size<Array<T, N>> : std::integral_constant<std::size_t, N> {};

template <std::size_t I, class T, std::size_t N>
auto get(Array<T, N> const &a) -> std::enable_if_t<(I < N), const T &> {
  return a[I];
}

} // namespace Utils
#endif
