#ifndef ESPRESSO_MASK_HPP
#define ESPRESSO_MASK_HPP

#include <utils/Array.hpp>
#include <utils/get.hpp>

#include <type_traits>
#include <utility>

namespace Utils {
namespace detail {
template <class T, class Integral, size_t... I>
auto mask_impl(Integral mask, T t, std::index_sequence<I...>) {
  return T{((mask & (1u << I)) ? get<I>(t) : tuple_element_t<I, T>{})...};
}
} // namespace detail

/**
 * @brief Set elements of a tuple-like to zero by a bit mask.
 *
 * @tparam T implements the tuple interface(get, tuple_size, ...)
 * @tparam Integral An unsigned integral type
 * @param mask bit mask, if the i-th bit is set, the i-th element
 *        in @param t is set to zero.
 * @param t
 * @return t partially zeroed out according to mask
 */
template <class T, class Integral>
std::enable_if_t<std::is_unsigned<Integral>::value &&
                     (8 * sizeof(Integral) >= tuple_size<T>::value),
                 T>
mask(Integral mask, T t) {
  return detail::mask_impl(mask, t, std::make_index_sequence<tuple_size<T>::value>{});
}
} // namespace Utils

#endif // ESPRESSO_MASK_HPP
