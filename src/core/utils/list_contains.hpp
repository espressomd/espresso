#ifndef UTILS_LIST_CONTAINS_HPP
#define UTILS_LIST_CONTAINS_HPP

#include "List.hpp"

namespace Utils {
/** @brief Check wether an @ref Utils::List contains the value c. */
template <typename T> bool list_contains(List<T> const &l, T const &c) {
  return std::any_of(l.begin(), l.end(), [c](T const &e) { return e == c; });
}
}

#endif
