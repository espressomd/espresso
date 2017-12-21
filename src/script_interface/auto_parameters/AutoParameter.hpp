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
/* Base case */
template <typename T> struct infer_length_helper {
  static constexpr size_t value{0};
};

/* Specialization for Vectors */
template <size_t N, typename T> struct infer_length_helper<Vector<N, T>> {
  static constexpr size_t value{N};
};
}

/**
 * @brief Infer supposed length of the paramer.
 *
 * This currently only works for fixed-sized vectors,
 * where the length is encoded in the type.
 */
template <typename T> constexpr size_t infer_length() {
  return detail::infer_length_helper<T>::value;
}

/**
 * @brief Describtion and getter/setter for a parameter.
 *
 * This is to be used with @c AutoParameters, see there for
 * more detailed documentation.
 */
struct AutoParameter {
  /* Exception types */
  struct WriteError {};

  /**
   * @brief read-write parameter that is bound to an object
   *        with setter and getter that are member functions.
   *
   * @param name The name the parameter should be bound to in the interface.
   * @param type The parameter type, by default this is deduced from the
   *             type of the reference.
   * @param length The supposed length of the parameter, by default this this
   *               is deduced from the type of the reference.
   */

  /* setter-getter */
  template <typename T, class O>
  AutoParameter(std::string const &name, std::shared_ptr<O> &obj,
                void (O::*setter)(T const &), T const &(O::*getter)() const,
                VariantType type = infer_type<T>(),
                size_t length = infer_length<T>())
      : name(name), type(type), length(length),
        set([&obj, setter](Variant const &v) {
          (obj.get()->*setter)(get_value<T>(v));
        }),
        get([&obj, getter]() { return (obj.get()->*getter)(); }) {}

  template <typename T, class O>
  AutoParameter(std::string const &name, std::shared_ptr<O> &obj,
                void (O::*setter)(T const &), T (O::*getter)() const,
                VariantType type = infer_type<T>(),
                size_t length = infer_length<T>())
      : name(name), type(type), length(length),
        set([&obj, setter](Variant const &v) {
          (obj.get()->*setter)(get_value<T>(v));
        }),
        get([&obj, getter]() { return (obj.get()->*getter)(); }) {}

  /* read-only */
  template <typename T, class O>
  AutoParameter(std::string const &name, std::shared_ptr<O> &obj,
                T const &(O::*getter)() const,
                VariantType type = infer_type<T>(),
                size_t length = infer_length<T>())
      : name(name), type(type), length(length),
        set([](Variant const &) { throw WriteError{}; }),
        get([&obj, getter]() { return (obj.get()->*getter)(); }) {}

  template <typename T, class O>
  AutoParameter(std::string const &name, std::shared_ptr<O> &obj,
                T (O::*getter)() const, VariantType type = infer_type<T>(),
                size_t length = infer_length<T>())
      : name(name), type(type), length(length),
        set([](Variant const &) { throw WriteError{}; }),
        get([&obj, getter]() { return (obj.get()->*getter)(); }) {}

  template <typename T, class O>
  AutoParameter(std::string const &name, std::shared_ptr<O> &obj, T O::*getter,
                VariantType type = infer_type<T>(),
                size_t length = infer_length<T>())
      : name(name), type(type), length(length),
        set([](Variant const &) { throw WriteError{}; }),
        get([&obj, getter]() { return (obj.get()->*getter)(); }) {}

  /* read-write */
  template <typename T, class O>
  AutoParameter(std::string const &name, std::shared_ptr<O> &obj,
                T &(O::*getter_setter)(), VariantType type = infer_type<T>(),
                size_t length = infer_length<T>())
      : name(name), type(type), length(length),
        set([&obj, getter_setter](Variant const &v) {
          (obj.get()->*getter_setter)() = get_value<T>(v);
        }),
        get([&obj, getter_setter]() { return (obj.get()->*getter_setter)(); }) {
  }

  /**
   * @brief read-write parameter that is bound to a referece.
   *
   * @param name The name the parameter should be bound to in the interface.
   * @param type The parameter type, by default this is deduced from the
   *             type of the reference.
   * @param length The supposed length of the parameter, by default this this
   *               is deduced from the type of the reference.
   */
  template <typename T>
  AutoParameter(std::string const &name, T &binding,
                VariantType type = infer_type<T>(),
                size_t length = infer_length<T>())
      : name(name), type(type), length(length),
        set([&binding](Variant const &v) { binding = get_value<T>(v); }),
        get([&binding]() { return Variant{binding}; }) {}

  template <typename T>
  AutoParameter(std::string const &name, std::shared_ptr<T> &binding,
                VariantType type = infer_type<std::shared_ptr<T>>(),
                size_t length = infer_length<std::shared_ptr<T>>())
      : name(name), type(type), length(length),
        set([&binding](Variant const &v) {
          binding = get_value<std::shared_ptr<T>>(v);
        }),
        get([&binding]() { return (binding) ? binding->id() : ObjectId(); }) {}

  /**
   * @brief read-only parameter that is bound to a const referece.
   *
   * @param name The name the parameter should be bound to in the interface.
   * @param type The parameter type, by default this is deduced from the
   *             type of the reference.
   * @param length The supposed length of the parameter, by default this this
   *               is deduced from the type of the reference.
   */
  template <typename T>
  AutoParameter(std::string const &name, T const &binding,
                VariantType type = infer_type<T>(),
                size_t length = infer_length<T>())
      : name(name), type(type), length(length),
        set([](Variant const &) { throw WriteError{}; }),
        get([&binding]() -> Variant { return binding; }) {}

  template <typename T>
  AutoParameter(std::string const &name, std::shared_ptr<T> const &binding,
                VariantType type = infer_type<std::shared_ptr<T>>(),
                size_t length = infer_length<std::shared_ptr<T>>())
      : name(name), type(type), length(length),
        set([](Variant const &) { throw WriteError{}; }),
        get([&binding]() { return (binding) ? binding->id() : ObjectId(); }) {}

  /**
   * @brief Parameter with a user-proivded getter and setter.
   *
   * @param name The name the parameter should be bound to in the interface.
   * @param set A setter, which can be a Functor, a Lambda or a std::function
   *            that accepts a @c Variant const&.
   * @param get A getter, which can be a Functor, a Lambda or a std::function
   *            that return the parameter. The return type must be convertible
   *            to Variant.
   * @param type The parameter type, by default this is deduced from the return
   *             type of the getter.
   * @param length The supposed length of the parameter, by default this this
   *               is deduced from the return type of the getter.
   */
  template <typename F, typename G,
            /* Try to guess the type from the return type of the getter */
            typename R = typename decltype(
                Utils::make_function(std::declval<G>()))::result_type>
  AutoParameter(std::string const &name, F const &set, G const &get,
                VariantType type = infer_type<R>(),
                size_t length = infer_length<R>())
      : name(name), type(type), length(length), set(Utils::make_function(set)),
        get(Utils::make_function(get)) {}

  /** The interface name. */
  const std::string name;
  /** The expected type. */
  VariantType type;
  /** The expected length. */
  size_t length;

  /**
   * @brief Set the parameter.
   */
  const std::function<void(Variant const &)> set;
  /**
   * @brief Get the parameter.
   */
  const std::function<Variant()> get;
};
}

#endif
