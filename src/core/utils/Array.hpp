#ifndef UTILS_ARRAY_HPP
#define UTILS_ARRAY_HPP

#include <cstddef>

#if defined(__CUDACC__)
#define DEVICE_QUALIFIER __host__ __device__
#else
#define DEVICE_QUALIFIER
#endif

namespace Utils {

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
  value_type m_data[N ? N : 1];
  DEVICE_QUALIFIER constexpr reference at(size_type i) {
    if (i >= N)
      throw std::out_of_range("Array access out of bounds.");
    return m_data[i];
  }
  DEVICE_QUALIFIER constexpr const_reference at(size_type i) const {
    if (i >= N)
      throw std::out_of_range("Array access out of bounds.");
    return m_data[i];
  }
  DEVICE_QUALIFIER constexpr reference operator[](size_type i) {
    return m_data[i];
  }
  DEVICE_QUALIFIER constexpr const_reference operator[](size_type i) const {
    return m_data[i];
  }
  DEVICE_QUALIFIER constexpr reference front() { return *begin(); }
  DEVICE_QUALIFIER constexpr const_reference front() const { return *cbegin(); }
  DEVICE_QUALIFIER constexpr reference back() {
    return N ? *(end() - 1) : *end();
  }
  DEVICE_QUALIFIER constexpr const_reference back() const {
    return N ? *(cend() - 1) : *cend();
  }
  DEVICE_QUALIFIER constexpr pointer data() noexcept { return &m_data[0]; }
  DEVICE_QUALIFIER constexpr const_pointer data() const noexcept {
    return &m_data[0];
  }
  DEVICE_QUALIFIER constexpr iterator begin() { return &m_data[0]; };
  DEVICE_QUALIFIER constexpr const_iterator cbegin() const {
    return &m_data[0];
  };
  DEVICE_QUALIFIER constexpr iterator end() { return &m_data[N]; };
  DEVICE_QUALIFIER constexpr const_iterator cend() const { return &m_data[N]; };
  DEVICE_QUALIFIER constexpr bool empty() const { return size() == 0; }
  DEVICE_QUALIFIER constexpr size_type size() const { return N; }
  DEVICE_QUALIFIER constexpr size_type max_size() const { return N; }
  DEVICE_QUALIFIER void fill(const value_type &value) {
    for (int i = 0; i < N; ++i)
      m_data[i] = value;
  }
};

} // namespace Utils
#undef DEVICE_QUALIFIER
#endif
