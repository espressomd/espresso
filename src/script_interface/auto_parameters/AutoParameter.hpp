/*
 * Copyright (C) 2010-2022 The ESPResSo project
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef SCRIPT_INTERFACE_AUTO_PARAMETERS_AUTO_PARAMETER_HPP
#define SCRIPT_INTERFACE_AUTO_PARAMETERS_AUTO_PARAMETER_HPP

#include "script_interface/Variant.hpp"
#include "script_interface/get_value.hpp"

#include <functional>
#include <memory>
#include <string>
#include <utility>

namespace ScriptInterface {
/**
 * @brief Description and getter/setter for a parameter.
 *
 * This is to be used with @c AutoParameters.
 */
struct AutoParameter {
  /* Exception types */
  struct WriteError {};

  /* Result types */
  struct ReadOnly {};
  static constexpr const ReadOnly read_only = ReadOnly{};

  /** @brief Read-write parameter that is bound to an object.
   *
   *  @param name    The name the parameter should be bound to in the interface.
   *  @param obj     The object the parameter should be bound to.
   *  @param setter  The setter.
   *  @param getter  The getter.
   */
  template <typename T, class O>
  AutoParameter(const char *name, std::shared_ptr<O> &obj,
                void (O::*setter)(T const &), T const &(O::*getter)() const)
      : name(name), setter_([&obj, setter](Variant const &v) {
          (obj.get()->*setter)(get_value<T>(v));
        }),
        getter_([&obj, getter]() { return (obj.get()->*getter)(); }) {}

  /** @brief Read-write parameter that is bound to an object.
   *  @overload
   */
  template <typename T, class O>
  AutoParameter(const char *name, std::shared_ptr<O> &obj,
                void (O::*setter)(T const &), T (O::*getter)() const)
      : name(name), setter_([&obj, setter](Variant const &v) {
          (obj.get()->*setter)(get_value<T>(v));
        }),
        getter_([&obj, getter]() { return (obj.get()->*getter)(); }) {}

  /** @brief Read-only parameter that is bound to an object.
   *
   *  @param name    The name the parameter should be bound to in the interface.
   *  @param obj     The object the parameter should be bound to.
   *  @param getter  The getter.
   */
  template <typename T, class O>
  AutoParameter(const char *name, std::shared_ptr<O> &obj,
                T const &(O::*getter)())
      : name(name), setter_([](Variant const &) { throw WriteError{}; }),
        getter_([&obj, getter]() { return (obj.get()->*getter)(); }) {}

  /** @brief Read-only parameter that is bound to an object.
   *  @overload
   */
  template <typename T, class O>
  AutoParameter(const char *name, std::shared_ptr<O> &obj,
                T (O::*getter)() const)
      : name(name), setter_([](Variant const &) { throw WriteError{}; }),
        getter_([&obj, getter]() { return (obj.get()->*getter)(); }) {}

  /** @brief Read-only parameter that is bound to an object.
   *  @overload
   */
  template <typename T, class O>
  AutoParameter(const char *name, std::shared_ptr<O> &obj, T O::*getter)
      : name(name), setter_([](Variant const &) { throw WriteError{}; }),
        getter_([&obj, getter]() { return (obj.get()->*getter)(); }) {}

  /** @brief Read-write parameter that is bound to an object.
   *  @overload
   */
  template <typename T, class O>
  AutoParameter(const char *name, std::shared_ptr<O> &obj,
                T &(O::*getter_setter)())
      : name(name), setter_([&obj, getter_setter](Variant const &v) {
          (obj.get()->*getter_setter)() = get_value<T>(v);
        }),
        getter_([&obj, getter_setter]() {
          return (obj.get()->*getter_setter)();
        }) {}

  /** @brief Read-write parameter that is bound to an object.
   *
   *  @param name   The name the parameter should be bound to in the interface.
   *  @param binding  The reference the parameter should be bound to.
   */
  template <typename T>
  AutoParameter(const char *name, T &binding)
      : name(name),
        setter_([&binding](Variant const &v) { binding = get_value<T>(v); }),
        getter_([&binding]() { return Variant{binding}; }) {}

  /** @brief Read-only parameter that is bound to an object.
   *  @overload
   */
  template <typename T>
  AutoParameter(const char *name, T const &binding)
      : name(name), setter_([](Variant const &) { throw WriteError{}; }),
        getter_([&binding]() -> Variant { return binding; }) {}

  /**
   * @brief Read-write parameter with a user-provided getter and setter.
   *
   * @param name The name the parameter should be bound to in the interface.
   * @param set A setter, which can be a Functor, a Lambda or a std::function
   *            that accepts a @c Variant const&.
   * @param get A getter, which can be a Functor, a Lambda or a std::function
   *            that return the parameter. The return type must be convertible
   *            to Variant.
   */
  template <typename Setter, typename Getter>
  AutoParameter(const char *name, Setter const &set, Getter const &get)
      : name(name), setter_(set), getter_(get) {}

  /**
   * @brief Read-only parameter with a user-provided getter.
   * @overload
   */
  template <typename Getter>
  AutoParameter(const char *name, ReadOnly, Getter const &get)
      : name(name), setter_([](Variant const &) { throw WriteError{}; }),
        getter_(get) {}

  /** The interface name. */
  const std::string name;

  /**
   * @brief Set the parameter.
   */
  const std::function<void(Variant const &)> setter_;
  /**
   * @brief Get the parameter.
   */
  const std::function<Variant()> getter_;

  void set(Variant const &v) const { setter_(v); }

  Variant get() const { return getter_(); }
};
} // namespace ScriptInterface

#endif
