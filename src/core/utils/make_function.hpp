#ifndef UTILS_MAKE_FUNCTION_HPP
#define UTILS_MAKE_FUNCTION_HPP

#include <functional>

#include "type_traits.hpp"

namespace Utils {

namespace detail {
template <class FPtr> struct function_traits;

template <class T, class C> struct function_traits<T(C::*)> { typedef T type; };
}

/* Identity for std::function */
template <typename T, typename... Args>
std::function<T(Args...)> make_function(std::function<T(Args...)> const &f) {
  return f;
}

/* For lambdas and functors */
template <
    typename F,
    typename Signature = typename Utils::function_remove_const<
        typename detail::function_traits<decltype(&F::operator())>::type>::type>
std::function<Signature> make_function(F const &f) {
  return f;
}

} /* namespace Utils */

#endif
