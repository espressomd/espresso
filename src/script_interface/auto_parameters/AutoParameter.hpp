/*
Copyright (C) 2010-2018 The ESPResSo project

This file is part of ESPResSo.

ESPResSo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ESPResSo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef SCRIPT_INTERFACE_AUTO_PARAMETERS_AUTO_PARAMETER_HPP
#define SCRIPT_INTERFACE_AUTO_PARAMETERS_AUTO_PARAMETER_HPP

#include <functional>
#include <memory>
#include <utility>

#include "Variant.hpp"
#include "get_value.hpp"

#include "utils/make_function.hpp"

namespace ScriptInterface {

namespace detail {
/* Base case */
template <typename T> struct infer_length_helper {
  static constexpr size_t value{0};
};

/* Specialization for Vectors */
template <size_t N, typename T> struct infer_length_helper<Vector<T, N>> {
  static constexpr size_t value{N};
};

template <size_t N, size_t M, typename T>
struct infer_length_helper<Vector<Vector<T, N>, M>> {
  static constexpr size_t value{N * M};
};
} // namespace detail

/**
 * @brief Infer supposed length of the parameter.
 *
 * This currently only works for fixed-sized vectors,
 * where the length is encoded in the type.
 */
template <typename T> constexpr size_t infer_length() {
  return detail::infer_length_helper<T>::value;
}

/**
 * @brief Description and getter/setter for a parameter.
 *
 * This is to be used with @c AutoParameters.
 */
struct AutoParameter {
  /* Exception types */
  struct WriteError {};

  /* Tag types */
  struct ReadOnly {};
  static constexpr const ReadOnly read_only = ReadOnly{};

  /** @brief Read-write parameter that is bound to an object.
   *
   *  @param name    The name the parameter should be bound to in the interface.
   *  @param obj     The object the parameter should be bound to.
   *  @param setter  The setter.
   *  @param getter  The getter.
   *  @param type    The parameter type, by default this is deduced from the
   *                 type of the reference.
   *  @param length  The supposed length of the parameter, by default this this
   *                 is deduced from the type of the reference.
   */
  template <typename T, class O>
  AutoParameter(std::string name, std::shared_ptr<O> &obj,
                void (O::*setter)(T const &), T const &(O::*getter)() const,
                VariantType type = infer_type<T>(),
                size_t length = infer_length<T>())
      : name(std::move(name)), type(type), length(length),
        set([&obj, setter](Variant const &v) {
          (obj.get()->*setter)(get_value<T>(v));
        }),
        get([&obj, getter]() { return (obj.get()->*getter)(); }) {}

  /** @brief Read-write parameter that is bound to an object.
   *  @overload
   */
  template <typename T, class O>
  AutoParameter(std::string name, std::shared_ptr<O> &obj,
                void (O::*setter)(T const &), T (O::*getter)() const,
                VariantType type = infer_type<T>(),
                size_t length = infer_length<T>())
      : name(std::move(name)), type(type), length(length),
        set([&obj, setter](Variant const &v) {
          (obj.get()->*setter)(get_value<T>(v));
        }),
        get([&obj, getter]() { return (obj.get()->*getter)(); }) {}

  /** @brief Read-only parameter that is bound to an object.
   *
   *  @param name    The name the parameter should be bound to in the interface.
   *  @param obj     The object the parameter should be bound to.
   *  @param getter  The getter.
   *  @param type    The parameter type, by default this is deduced from the
   *                 type of the reference.
   *  @param length  The supposed length of the parameter, by default this this
   *                 is deduced from the type of the reference.
   */
  template <typename T, class O>
  AutoParameter(std::string name, std::shared_ptr<O> &obj,
                T const &(O::*getter)() const,
                VariantType type = infer_type<T>(),
                size_t length = infer_length<T>())
      : name(std::move(name)), type(type), length(length),
        set([](Variant const &) { throw WriteError{}; }),
        get([&obj, getter]() { return (obj.get()->*getter)(); }) {}

  /** @brief Read-only parameter that is bound to an object.
   *  @overload
   */
  template <typename T, class O>
  AutoParameter(std::string name, std::shared_ptr<O> &obj,
                T (O::*getter)() const, VariantType type = infer_type<T>(),
                size_t length = infer_length<T>())
      : name(std::move(name)), type(type), length(length),
        set([](Variant const &) { throw WriteError{}; }),
        get([&obj, getter]() { return (obj.get()->*getter)(); }) {}

  /** @brief Read-only parameter that is bound to an object.
   *  @overload
   */
  template <typename T, class O>
  AutoParameter(std::string name, std::shared_ptr<O> &obj, T O::*getter,
                VariantType type = infer_type<T>(),
                size_t length = infer_length<T>())
      : name(std::move(name)), type(type), length(length),
        set([](Variant const &) { throw WriteError{}; }),
        get([&obj, getter]() { return (obj.get()->*getter)(); }) {}

  /** @brief Read-write parameter that is bound to an object.
   *  @overload
   */
  template <typename T, class O>
  AutoParameter(std::string name, std::shared_ptr<O> &obj,
                T &(O::*getter_setter)(), VariantType type = infer_type<T>(),
                size_t length = infer_length<T>())
      : name(std::move(name)), type(type), length(length),
        set([&obj, getter_setter](Variant const &v) {
          (obj.get()->*getter_setter)() = get_value<T>(v);
        }),
        get([&obj, getter_setter]() { return (obj.get()->*getter_setter)(); }) {
  }

  /** @brief Read-write parameter that is bound to an object.
   *
   *  @param name   The name the parameter should be bound to in the interface.
   *  @param binding  The reference the parameter should be bound to.
   *  @param type   The parameter type, by default this is deduced from the
   *                type of the reference.
   *  @param length The supposed length of the parameter, by default it
   *                is deduced from the type of the reference.
   */
  template <typename T>
  AutoParameter(std::string name, T &binding,
                VariantType type = infer_type<T>(),
                size_t length = infer_length<T>())
      : name(std::move(name)), type(type), length(length),
        set([&binding](Variant const &v) { binding = get_value<T>(v); }),
        get([&binding]() { return Variant{binding}; }) {}

  /** @brief Read-write parameter that is bound to an object.
   *  @overload
   */
  template <typename T>
  AutoParameter(std::string name, std::shared_ptr<T> &binding,
                VariantType type = infer_type<std::shared_ptr<T>>(),
                size_t length = infer_length<std::shared_ptr<T>>())
      : name(std::move(name)), type(type), length(length),
        set([&binding](Variant const &v) {
          binding = get_value<std::shared_ptr<T>>(v);
        }),
        get([&binding]() { return (binding) ? binding->id() : ObjectId(); }) {}

  /** @brief Read-only parameter that is bound to an object.
   *  @overload
   */
  template <typename T>
  AutoParameter(std::string name, T const &binding,
                VariantType type = infer_type<T>(),
                size_t length = infer_length<T>())
      : name(std::move(name)), type(type), length(length),
        set([](Variant const &) { throw WriteError{}; }),
        get([&binding]() -> Variant { return binding; }) {}

  /** @brief Read-only parameter that is bound to an object.
   *  @overload
   */
  template <typename T>
  AutoParameter(std::string name, std::shared_ptr<T> const &binding,
                VariantType type = infer_type<std::shared_ptr<T>>(),
                size_t length = infer_length<std::shared_ptr<T>>())
      : name(std::move(name)), type(type), length(length),
        set([](Variant const &) { throw WriteError{}; }),
        get([&binding]() { return (binding) ? binding->id() : ObjectId(); }) {}

  /**
   * @brief Read-write parameter with a user-provided getter and setter.
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
  AutoParameter(std::string name, F const &set, G const &get,
                VariantType type = infer_type<R>(),
                size_t length = infer_length<R>())
      : name(std::move(name)), type(type), length(length),
        set(Utils::make_function(set)), get(Utils::make_function(get)) {}

  /**
   * @brief Read-only parameter with a user-provided getter.
   * @overload
   */
  template <typename G,
            /* Try to guess the type from the return type of the getter */
            typename R = typename decltype(
                Utils::make_function(std::declval<G>()))::result_type>
  AutoParameter(std::string name, ReadOnly, G const &get,
                VariantType type = infer_type<R>(),
                size_t length = infer_length<R>())
      : name(std::move(name)), type(type), length(length),
        set([](Variant const &) { throw WriteError{}; }),
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
} // namespace ScriptInterface

#endif
