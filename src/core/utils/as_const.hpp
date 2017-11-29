#ifndef CORE_UTILS_AS_CONST_HPP
#define CORE_UTILS_AS_CONST_HPP

#include <type_traits>

namespace Utils {
template <class T>
constexpr typename std::add_const<T>::type &as_const(T &t) noexcept {
  return t;
}
}

#endif
