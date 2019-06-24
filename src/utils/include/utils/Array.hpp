#ifndef UTILS_ARRAY_HPP
#define UTILS_ARRAY_HPP

#include "device_qualifier.hpp"

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
    if (i >= N)
      DEVICE_THROW(std::out_of_range("Array access out of bounds."));
    return m_storage.m_data[i];
  }

  DEVICE_QUALIFIER constexpr const_reference at(size_type i) const {
    if (i >= N)
      DEVICE_THROW(std::out_of_range("Array access out of bounds."));
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
  __attribute__((annotate("serialize"))) void
  __cxxamp_serialize(Kalmar::Serialize &s) const {
    for (int i = 0; i < N; ++i)
      s.Append(sizeof(T), &m_storage[i]);
    for (int i = N; i < 10; ++i)
      s.Append(sizeof(T), T());
  }

  __attribute__((annotate("user_deserialize"))) void
  Foo(T x0, T x1, T x2, T x3, T x4, T x5, T x6, T x7, T x8, T x9)[[cpu]][[hc]] {
    if (N > 0)
      m_storage[0] = x0;
    if (N > 1)
      m_storage[1] = x1;
    if (N > 2)
      m_storage[2] = x2;
    if (N > 3)
      m_storage[3] = x3;
    if (N > 4)
      m_storage[4] = x4;
    if (N > 5)
      m_storage[5] = x5;
    if (N > 6)
      m_storage[6] = x6;
    if (N > 7)
      m_storage[7] = x7;
    if (N > 8)
      m_storage[8] = x8;
    if (N > 9)
      m_storage[9] = x9;
    static_assert(N <= 10, "custom serialization function is only implemented "
                           "for arrays of sizes up to 10");
  }
#endif

private:
  friend boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, const unsigned int /* version */) {
    ar &m_storage;
  }
};
} // namespace Utils
#endif
