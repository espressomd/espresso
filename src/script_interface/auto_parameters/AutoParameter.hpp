#ifndef SCRIPT_INTERFACE_AUTO_PARAMETERS_AUTO_PARAMETER_HPP
#define SCRIPT_INTERFACE_AUTO_PARAMETERS_AUTO_PARAMETER_HPP

#include <functional>
#include <memory>

#include "script_interface/Variant.hpp"

namespace ScriptInterface {
struct AutoParameter {
  /* read-write parameter */
  template <typename T>
  AutoParameter(std::string const &name, T &binding)
      : name(name), type(infer_type<T>()),
        set([&binding](Variant const &v) { binding = boost::get<T>(v); }),
        get([&binding]() -> Variant { return binding; }) {}

  /* read-only parameter */
  template <typename T>
  AutoParameter(std::string const &name, T const &binding)
      : name(name), type(infer_type<T>()),
        get([&binding]() -> Variant { return binding; }) {}

  /* user-provided */
  AutoParameter(std::string const &name, VariantType type,
                std::function<void(Variant const &)> set,
                std::function<Variant()> get)
      : name(name), type(type), set(set), get(get) {}

  const std::string name;
  VariantType type;
  const std::function<void(Variant const &)> set;
  const std::function<Variant()> get;
};
}

#endif
