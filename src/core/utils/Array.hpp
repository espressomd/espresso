#ifndef UTILS_ARRAY_HPP
#define UTILS_ARRAY_HPP

#include "device_qualifier.hpp"

#include <boost/serialization/access.hpp>
#include <cstddef>

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
      throw std::out_of_range("Array access out of bounds.");
    return m_storage.m_data[i];
  }

  DEVICE_QUALIFIER constexpr const_reference at(size_type i) const {
    if (i >= N)
      throw std::out_of_range("Array access out of bounds.");
    return m_storage.m_data[i];
  }

  DEVICE_QUALIFIER constexpr reference operator[](size_type i) {
    assert(i < N);
    return m_storage.m_data[i];
  }

  DEVICE_QUALIFIER constexpr const_reference operator[](size_type i) const {
    assert(i < N);
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

private:
  friend boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, const unsigned int /* version */) {
    ar &m_storage;
  }
};
} // namespace Utils
#endif
