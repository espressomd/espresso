#ifndef UTILS_TYPE_TRAITS_HPP
#define UTILS_TYPE_TRAITS_HPP

#include <type_traits>

namespace Utils {

/**
 * @brief Remove const from a function signature.
 */
template <typename T> struct function_remove_const;

template <typename R, typename... Args>
struct function_remove_const<R(Args...)> {
  using type = R(Args...);
};

template <typename R, typename... Args>
struct function_remove_const<R(Args...) const> {
  using type = R(Args...);
};

/**
 * @brief True iff T is an innstantiation of of Template.
 */
template <typename T, template <typename...> class Template>
struct is_instance_of : public std::false_type {};

template <typename... T, template <typename...> class Template>
struct is_instance_of<Template<T...>, Template> : public std::true_type {};
} // namespace Utils

#endif
