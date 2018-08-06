#ifndef UTILS_SPAN_HPP
#define UTILS_SPAN_HPP

#include <cassert>
#include <cstddef>
#include <iterator>
#include <stdexcept>
#include <type_traits>

namespace Utils {
/**
 * @brief A sripped-down version of std::span from C++17.
 *
 * Behaves like a std::span where implemented
 */

template <class T> class Span {
public:
  using value_type = typename std::remove_cv<T>::type;
  using pointer = T *;
  using const_pointer = const T *;
  using reference = T &;
  using const_reference = const T &;
  using iterator = pointer;
  using const_iterator = const_pointer;
  using reverse_iterator = std::reverse_iterator<iterator>;
  using const_reverse_iterator = std::reverse_iterator<const_iterator>;
  using size_type = size_t;
  using difference_type = ptrdiff_t;

private:
  T *m_ptr;
  size_t m_size;

public:
  Span() = default;
  Span(const Span &) = default;
  Span &operator=(const Span &) = default;

  Span(pointer array, size_type length) : m_ptr(array), m_size(length) {}

  size_type size() const { return m_size; }
  bool empty() const { return size() == 0; }

  iterator begin() const { return m_ptr; }
  const_iterator cbegin() const { return m_ptr; }
  iterator end() const { return m_ptr + m_size; }
  const_iterator cend() const { return m_ptr + m_size; }

  reference operator[](size_type i) const {
    assert(i < size());
    return m_ptr[i];
  }

  reference at(size_type i) const {
    return (i < size())
               ? m_ptr[i]
      : throw std::out_of_range("span access out of bounds."), m_ptr[i];
  }

  pointer data() const { return m_ptr; }
};
}

#endif
