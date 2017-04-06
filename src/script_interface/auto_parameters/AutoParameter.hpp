#ifndef SCRIPT_INTERFACE_AUTO_PARAMETERS_AUTO_PARAMETER_HPP
#define SCRIPT_INTERFACE_AUTO_PARAMETERS_AUTO_PARAMETER_HPP

#include <functional>
#include <memory>
#include <utility>

#include "core/utils/make_function.hpp"
#include "script_interface/Variant.hpp"
#include "script_interface/get_value.hpp"

namespace ScriptInterface {

namespace detail {
template <typename T> struct infer_length_helper {
  static constexpr size_t value{0};
};

template <size_t N, typename T> struct infer_length_helper<Vector<N, T>> {
  static constexpr size_t value{N};
};
}

template <typename T> constexpr size_t infer_length() {
  return detail::infer_length_helper<T>::value;
}

struct AutoParameter {
  /* Exception types */
  struct WriteError {};

  /* read-write parameter */
  template <typename T>
  AutoParameter(std::string const &name, T &binding,
                VariantType type = infer_type<T>(),
                size_t length = infer_length<T>())
      : name(name), type(type), length(length),
        set([&binding](Variant const &v) { binding = get_value<T>(v); }),
        get([&binding]() { return binding; }) {}

  /* read-only parameter */
  template <typename T>
  AutoParameter(std::string const &name, T const &binding,
                VariantType type = infer_type<T>(),
                size_t length = infer_length<T>())
      : name(name), type(type), length(length),
        set([this](Variant const &) { throw WriteError{}; }),
        get([&binding]() { return binding; }) {}

  /* user-provided getter and setter */
  template <typename F, typename G,
            /* Try to guess the type from the return type of the getter */
            typename R = typename decltype(
                Utils::make_function(std::declval<G>()))::result_type>
  AutoParameter(std::string const &name, F const &set, G const &get,
                VariantType type = infer_type<R>(),
                size_t length = infer_length<R>())
      : name(name), type(type), length(length), set(Utils::make_function(set)),
        get(Utils::make_function(get)) {}

  const std::function<void(Variant const &)> set;
  const std::function<Variant()> get;

  const std::string name;
  VariantType type;
  size_t length;
};
}

#endif
