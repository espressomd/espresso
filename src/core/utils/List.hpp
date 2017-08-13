#ifndef CORE_UTILS_LIST_HPP
#define CORE_UTILS_LIST_HPP

#include <cassert>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <iterator>
#include <algorithm>

#include "memory.hpp"

namespace Utils {
/** List.
    Use the functions specified in list operations. */
template <typename T> struct List {
  static_assert(std::is_pod<T>::value, "");

  using size_type = uint32_t;

  List() : e{nullptr}, n{0}, max{0} {}
  explicit List(size_type size) : List() { resize(n = size); }
  ~List() { resize(0); }

  template <size_t N> explicit List(T (&array)[N]) : List(N) {
    std::copy(std::begin(array), std::end(array), begin());
  }

private:
  void copy(List const &rhs) {
    resize(rhs.n);
    std::copy(rhs.begin(), rhs.end(), begin());
  }

  void move(List &&rhs) {
    using std::swap;
    swap(n, rhs.n);
    swap(max, rhs.max);
    swap(e, rhs.e);
  }

public:
  List(List const &rhs) : List() { copy(rhs); }
  List(List &&rhs) : List() { move(std::move(rhs)); }
  List &operator=(List const &rhs) {
    copy(rhs);
    return *this;
  }
  List &operator=(List &&rhs) {
    move(std::move(rhs));
    return *this;
  }

  T *begin() { return e; }
  T const *begin() const { return e; }
  T *end() { return e + n; }
  T const *end() const { return e + n; }
  size_type size() const { return n; }
  bool empty() const { return n == 0; }
  size_type capacity() const { return max; }
  void clear() { resize(0); }

private:
  /**
   * @brief Realloc memory in an exception safe way.
   *
   * If Utils::realloc fails, the original memory block
   * is unchanged an still vaild, but Utils::realloc will
   * throw. Because this->e is then not updated the List
   * actually stays unchanged, so that
   * we can give the strong exception safety guarantee here.
   */
  void realloc(size_type size) {
    auto new_ptr = Utils::realloc(this->e, sizeof(T) * size);
    this->e = new_ptr;
    this->max = size;
  }

public:
  void reserve(size_type size) {
    assert(size <= max_size());
    if (size < this->max) {
      realloc(size);
    }
  }

  /**
   * @brief Resize the list
   *
   * If the size is smaller than the current
   * size, the excess memory is free'd. This
   * is different from std::vector. If size
   * is bigger than the capacity, the new memory
   * is uninitialized.
   */
  void resize(size_type size) {
    assert(size <= max_size());
    if (size != this->max) {
      realloc(size);
      this->n = size;
    }
  }

public:
  size_type max_size() const { return std::numeric_limits<size_type>::max(); }

  T &operator[](size_type i) {
    assert(i < n);
    return e[i];
  }

  T const &operator[](size_type i) const {
    assert(i < n);
    return e[i];
  }

  /** Dynamically allocated field. */
  T *e;
  /** number of used elements in the field. */
  size_type n;
  /** allocated size of the field. This value is ONLY changed
      in the routines specified in list operations ! */
  size_type max;
};

} /* namespace Utils */

using IntList = Utils::List<int>;
using DoubleList = Utils::List<double>;

#endif
