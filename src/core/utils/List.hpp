#ifndef CORE_UTILS_LIST_HPP
#define CORE_UTILS_LIST_HPP

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iterator>
#include <limits>
#include <type_traits>

#include "memory.hpp"

namespace Utils {
/** List.
    Use the functions specified in list operations. */
template <typename T, typename SizeType = uint32_t> class List {
  static_assert(std::is_unsigned<SizeType>::value, "SizeType needs to be unsigned.");

public:
  using value_type = T;
  using size_type = SizeType;
  using difference_type = typename std::make_signed<size_type>::type;
  using iterator = T *;
  using const_iterator = T const *;
  using pointer = T *;
  using reference = T &;

public:
  List() noexcept : e{nullptr}, n{0}, max{0} {}
  explicit List(size_type size) : List() { resize(size); }
  List(size_type size, T const &value) : List(size) {
    std::fill(begin(), end(), value);
  }
  ~List() { resize(0); }

  template <size_t N> explicit List(T (&array)[N]) : List(N) {
    std::copy(std::begin(array), std::end(array), begin());
  }

  List(std::initializer_list<T> il) : List(il.size()) {
    std::copy(std::begin(il), std::end(il), begin());
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
  List(List &&rhs) noexcept : List() { move(std::move(rhs)); }
  List &operator=(List const &rhs) {
    copy(rhs);
    return *this; // NOLINT
  }
  List &operator=(List &&rhs) noexcept {
    move(std::move(rhs));
    return *this; // NOLINT
  }

  T *begin() { return e; }
  T const *begin() const { return e; }
  T *end() { return e + n; }
  T const *end() const { return e + n; }

  T &front() { return this->e[0]; }
  T const &front() const { return this->e[0]; }

  T &back() { return this->e[size() - 1]; }
  T const &back() const { return this->e[size() - 1]; }

  size_type size() const { return n; }
  bool empty() const { return n == 0; }
  size_type capacity() const { return max; }
  void clear() { resize(0); }
  T *data() { return e; }
  T const *data() const { return e; }

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
    if (size > this->max) {
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
    if (size != capacity()) {
      realloc(size);
      this->n = size;
    }
  }

  void push_back(T const &v) {
    auto const new_size = size() + 1;

    if (new_size > capacity()) {
      resize(new_size);
    } else {
      this->n = new_size;
    }

    this->back() = v;
  }

  template <typename... Args> void emplace_back(Args &&... args) {
    auto const new_size = size() + 1;

    if (new_size > capacity()) {
      resize(new_size);
    } else {
      this->n = new_size;
    }

    new (&(this->back())) T(std::forward<Args>(args)...);
  }

  /**
   * @brief Erase elements [first, last).
   */
  iterator erase(iterator first, iterator last) {
    auto const n_elem = std::distance(first, last);
    assert(n_elem >= 0);

    auto r = std::copy(last, end(), first);

    this->n -= n_elem;
    return --r;
  }

  void shrink_to_fit() {
    if (size() < capacity()) {
      resize(size());
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

  bool operator==(List const &rhs) const {
    return (size() == rhs.size()) &&
           std::equal(rhs.begin(), rhs.end(), begin());
  }

  bool operator!=(List const &rhs) const { return not(this->operator==(rhs)); }

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
