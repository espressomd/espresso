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
#ifndef SRC_UTILS_INCLUDE_UTILS_ARRAY_HPP
#define SRC_UTILS_INCLUDE_UTILS_ARRAY_HPP

/**
 * @file
 *
 * @brief Array implementation with CUDA support.
 */

#include "device_qualifier.hpp"
#include "get.hpp"
#include "serialization/array.hpp"

#include <boost/serialization/access.hpp>

#include <cassert>
#include <cstddef>
#include <iterator>
#include <ostream>
#include <stdexcept>

namespace Utils {
namespace detail {

template <typename T, std::size_t N> struct Storage {
  T m_data[N];

private:
  friend boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, const unsigned int /* version */) {
    ar & m_data;
  }
};

template <typename T> struct Storage<T, 0> {};

struct ArrayFormatterStream {
  std::ostream &stream;
  char const *separator;
  ArrayFormatterStream(std::ostream &s, char const *sep)
      : stream(s), separator(sep) {}
};

struct ArrayFormatter {
  char const *separator;
  friend ArrayFormatterStream operator<<(std::ostream &os,
                                         ArrayFormatter const &fmt) {
    return {os, fmt.separator};
  }
};

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
  }

  DEVICE_QUALIFIER constexpr const_iterator begin() const noexcept {
    return &m_storage.m_data[0];
  }

  DEVICE_QUALIFIER constexpr const_iterator cbegin() const noexcept {
    return &m_storage.m_data[0];
  }

  DEVICE_QUALIFIER constexpr iterator end() noexcept {
    return &m_storage.m_data[N];
  }

  DEVICE_QUALIFIER constexpr const_iterator end() const noexcept {
    return &m_storage.m_data[N];
  }

  DEVICE_QUALIFIER constexpr const_iterator cend() const noexcept {
    return &m_storage.m_data[N];
  }

  DEVICE_QUALIFIER constexpr bool empty() const noexcept { return size() == 0; }

  DEVICE_QUALIFIER constexpr size_type size() const noexcept { return N; }

  DEVICE_QUALIFIER constexpr size_type max_size() const noexcept { return N; }

  DEVICE_QUALIFIER void fill(const value_type &value) {
    for (size_type i = 0; i != size(); ++i) {
      m_storage.m_data[i] = value;
    }
  }

  DEVICE_QUALIFIER static constexpr Array<T, N>
  broadcast(const value_type &value) {
    Array<T, N> ret{};
    for (size_type i = 0; i != N; ++i) {
      ret[i] = value;
    }
    return ret;
  }

  static constexpr detail::ArrayFormatter formatter(char const *sep = " ") {
    return {sep};
  }

private:
  friend boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, const unsigned int /* version */) {
    ar & m_storage;
  }

  static std::ostream &format(std::ostream &out, Array const &a,
                              char const *const sep) {
    if (not a.empty())
      out << a[0];

    for (auto it = std::next(a.begin()); it != a.end(); ++it)
      out << sep << *it;

    return out;
  }

  friend std::ostream &operator<<(std::ostream &out, Array const &a) {
    return format(out, a, ", ");
  }

  friend std::ostream &operator<<(detail::ArrayFormatterStream const &fmt,
                                  Array const &a) {
    return format(fmt.stream, a, fmt.separator);
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

UTILS_ARRAY_BOOST_MPI_T(Utils::detail::Storage, N)
UTILS_ARRAY_BOOST_BIT_S(Utils::detail::Storage, N)
UTILS_ARRAY_BOOST_CLASS(Utils::detail::Storage, N, object_serializable)
UTILS_ARRAY_BOOST_TRACK(Utils::detail::Storage, N, track_never)
UTILS_ARRAY_BOOST_MPI_T(Utils::Array, N)
UTILS_ARRAY_BOOST_BIT_S(Utils::Array, N)
UTILS_ARRAY_BOOST_CLASS(Utils::Array, N, object_serializable)
UTILS_ARRAY_BOOST_TRACK(Utils::Array, N, track_never)

#endif
