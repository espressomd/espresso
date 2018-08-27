#ifndef UTILS_MAKE_FUNCTION_HPP
#define UTILS_MAKE_FUNCTION_HPP

#include <functional>

#include "type_traits.hpp"

namespace Utils {

namespace detail {
template <class FPtr> struct function_traits;

/**
 * @brief Determine the signature from a pointer to member function (C::*).
 */
template <class T, class C> struct function_traits<T(C::*)> { typedef T type; };
} // namespace detail

/* make_function deduces the signature of a class with not-overloaded operator()
   member , a lambda or std::function,
   and creates a corresponding std::function. This is needed to solve some
   issues with
   template parameer deduction if the type is not now beforehand. */

/**
 * @brief Given a std::function, return a copy.
 *
 * For a std::function argument nothing has to be done, just a copy is returned.
 */
template <typename T, typename... Args>
std::function<T(Args...)> make_function(std::function<T(Args...)> const &f) {
  return f;
}

/**
 * @brief Given a Functor or a Closure, return a std::function with the same
 * signature that contains a copy of the argument.
 *
 * This inspects the operator() on the given type, deduces the signature,
 * if operator() is const, as in non-mutable lambdas, the const is removed,
 * and a std::function with the correct type, containing a copy, is returned.
 */
template <
    typename F,
    typename Signature = typename Utils::function_remove_const<
        typename detail::function_traits<decltype(&F::operator())>::type>::type>
std::function<Signature> make_function(F const &f) {
  return f;
}

} /* namespace Utils */

#endif
