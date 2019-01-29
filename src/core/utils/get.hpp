#ifndef UTILS_GET_HPP
#define UTILS_GET_HPP

namespace Utils {

template<std::size_t I, typename T>
const auto get(const T& v) {
  return std::get<I>(v);
}

}
#endif
