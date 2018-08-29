#ifndef UTILS_SPAN_HPP
#define UTILS_SPAN_HPP

#include <cassert>
#include <cstddef>
#include <iterator>
#include <stdexcept>
#include <type_traits>

namespace Utils {
namespace detail {
template <typename T> using decay_t = typename std::decay<T>::type;

template <class T, class C>
using has_data =
    std::is_convertible<decay_t<decltype(std::declval<C>().data())> *,
                        T *const *>;
} // namespace detail

/**
 * @brief A stripped-down version of std::span from C++17.
 *
 * Behaves like a std::span where implemented.
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
  size_t m_size{};

  template <typename U>
  using enable_if_const_t =
      typename std::enable_if<std::is_const<T>::value, U>::type;
  template <class U>
  using enable_if_mutable_t =
      typename std::enable_if<!std::is_const<T>::value, U>::type;
  template <class U>
  using enable_if_has_data_t =
      typename std::enable_if<detail::has_data<T, U>::value, U>::type;

public:
  Span() = default;
  Span(const Span &) = default;
  Span &operator=(const Span &) = default;

  constexpr Span(pointer array, size_type length)
      : m_ptr(array), m_size(length) {}
  template <size_t N> constexpr Span(T (&a)[N]) noexcept : Span(a, N) {}

  template <typename C, typename = enable_if_mutable_t<C>,
            typename = enable_if_has_data_t<C>>
  explicit Span(C &c) noexcept : Span(c.data(), c.size()) {}
  template <typename C, typename = enable_if_const_t<C>,
            typename = enable_if_has_data_t<C>>
  Span(const C &c) noexcept : Span(c.data(), c.size()) {}

  constexpr size_type size() const { return m_size; }
  constexpr bool empty() const { return size() == 0; }

  constexpr iterator begin() const { return m_ptr; }
  constexpr const_iterator cbegin() const { return m_ptr; }
  constexpr iterator end() const { return m_ptr + m_size; }
  constexpr const_iterator cend() const { return m_ptr + m_size; }
  constexpr reverse_iterator rbegin() const { return reverse_iterator(end()); }
  constexpr reverse_iterator rend() const { return reverse_iterator(begin()); }

  constexpr reference operator[](size_type i) const {
    return assert(i < size()), m_ptr[i];
  }

  constexpr reference at(size_type i) const {
    return (i < size()) ? m_ptr[i]
                        : throw std::out_of_range("span access out of bounds."),
           m_ptr[i];
  }

  constexpr pointer data() const { return m_ptr; }
};
} // namespace Utils

#endif
