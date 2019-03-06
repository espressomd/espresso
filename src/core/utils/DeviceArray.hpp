#ifndef UTILS_DEVICE_ARRAY_HPP
#define UTILS_DEVICE_ARRAY_HPP

#include <cstddef>

#if defined(__CUDACC__)
#define QUALIFIERS __host__ __device__
#else
#define QUALIFIERS
#endif 

namespace Utils {

template<typename T, std::size_t N>
struct DeviceArray {
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
  QUALIFIERS constexpr reference at(size_type i) {
    return (i < N) ? m_data[i] : throw std::out_of_range("DeviceArray access out of bounds.");
  }
  QUALIFIERS constexpr const_reference at(size_type i) const {
    return (i < N) ? m_data[i] : throw std::out_of_range("DeviceArray access out of bounds.");
  }
  QUALIFIERS constexpr reference operator[](size_type i) {
    return m_data[i];
  }
  QUALIFIERS constexpr const_reference operator[](size_type i) const {
    return m_data[i];
  }
  QUALIFIERS constexpr reference front() {
    return *begin();
  }
  QUALIFIERS constexpr const_reference front() const {
    return *cbegin();
  }
  QUALIFIERS constexpr reference back() {
    return N ? *(end() - 1) : *end();
  }
  QUALIFIERS constexpr const_reference back() const {
    return N ? *(cend() - 1) : *cend();
  }
  QUALIFIERS constexpr pointer data() noexcept {
    return &m_data[0];
  }
  QUALIFIERS constexpr const_pointer data() const noexcept {
    return &m_data[0];
  }
  QUALIFIERS constexpr iterator begin() { return &m_data[0]; };
  QUALIFIERS constexpr const_iterator cbegin() const { return &m_data[0]; };
  QUALIFIERS constexpr iterator end() { return &m_data[N]; };
  QUALIFIERS constexpr const_iterator cend() const { return &m_data[N]; };
  QUALIFIERS constexpr bool empty() const { return size() == 0; }
  QUALIFIERS constexpr size_type size() const { return N; }
  QUALIFIERS constexpr size_type max_size() const { return N; }
  QUALIFIERS void fill(const value_type& value) {
    for (int i=0; i<N; ++i) m_data[i] = value;
  }
};

} // namespace Utils
#undef QUALIFIERS
#endif
