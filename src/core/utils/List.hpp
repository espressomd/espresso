#ifndef CORE_UTILS_LIST_HPP
#define CORE_UTILS_LIST_HPP

#include "memory.hpp"

namespace Utils {

/** List.
    Use the functions specified in list operations. */
template <typename T> struct List {
  T *begin() { return e; }
  T const *begin() const { return e; }
  T *end() { return e + n; }
  T const *end() const { return e + n; }
  int size() const { return n; }
  void resize(int size) {
    if (size != this->max) {
      this->max = size;
      this->e = Utils::realloc(this->e, sizeof(int) * this->max);
    }
  }
  /** Dynamically allocated field. */
  T *e;
  /** number of used elements in the field. */
  int n;
  /** allocated size of the field. This value is ONLY changed
      in the routines specified in list operations ! */
  int max;
};

} /* namespace Utils */

using IntList = Utils::List<int>;
using DoubleList = Utils::List<double>;

#endif
