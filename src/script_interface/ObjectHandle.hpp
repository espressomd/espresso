/*
 * Copyright (C) 2015-2022 The ESPResSo project
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

#ifndef SCRIPT_INTERFACE_SCRIPT_INTERFACE_BASE_HPP
#define SCRIPT_INTERFACE_SCRIPT_INTERFACE_BASE_HPP
#include "Variant.hpp"

#include <utils/Span.hpp>

#include <boost/utility/string_ref.hpp>

#include <memory>
#include <string>
#include <vector>

namespace ScriptInterface {
class Context;

/**
 * @brief Base class for interface handles.
 */
class ObjectHandle {
public:
  ObjectHandle() = default;

  /* Copy has unclear semantics, so it should not be allowed. */
  ObjectHandle(ObjectHandle const &) = delete;
  ObjectHandle &operator=(ObjectHandle const &) = delete;
  virtual ~ObjectHandle() = default;

private:
  friend class Context;
  std::shared_ptr<Context> m_context = {};

public:
  boost::string_ref name() const;

public:
  /**
   * @brief Responsible context.
   */
  Context *context() const { return m_context.get(); }

public:
  /**
   * @brief Construct the handled object.
   *
   * This function is called on object creation with user
   * provided parameters. This can be used if the SO has required parameters,
   * it represents some type that can not reasonably be default constructed,
   * or if the core implementation has to be chosen by a parameter. It is
   * guaranteed that no getter or setter functions from this interface is
   * called before construction (only @ref name() and @ref valid_parameters()),
   * and it is only called once.
   *
   * The default implementation just calls @ref set_parameter for every
   * parameter.
   *
   * @param params The parameters to the constructor. Only parameters that
   *               are valid for a default-constructed object are valid.
   */
  void construct(VariantMap const &params) { do_construct(params); }

private:
  virtual void do_construct(VariantMap const &params) {
    for (auto const &p : params) {
      do_set_parameter(p.first, p.second);
    }
  }

public:
  /**
   * @brief Get current parameters.
   * @return Parameters set in class.
   */
  VariantMap get_parameters() const {
    VariantMap values;

    for (auto const &p : valid_parameters()) {
      values[p.data()] = get_parameter(p.data());
    }

    return values;
  }

  /**
   * @brief Get required and optional parameters for class.
   * @return Expected parameters.
   */
  virtual Utils::Span<const boost::string_ref> valid_parameters() const {
    return {};
  }

  auto get_valid_parameters() const {
    auto const names = valid_parameters();
    return std::vector<std::string>(names.begin(), names.end());
  }

  /**
   * @brief Get single parameter.
   *
   * @param name Name of the parameter
   * @return Value of parameter @p name
   */
  virtual Variant get_parameter(const std::string &name) const { return {}; }

  /**
   * @brief Set single parameter.
   */
  void set_parameter(const std::string &name, const Variant &value);

private:
  /**
   * @brief Local implementation of @ref set_parameter.
   */
  virtual void do_set_parameter(const std::string &, const Variant &) {}

public:
  /**
   * @brief Call a method on the object.
   */
  Variant call_method(const std::string &name, const VariantMap &params);

protected:
  /**
   * @brief Local implementation of @c call_method.
   *
   * If not overridden by the implementation, this does nothing.
   */
  virtual Variant do_call_method(const std::string &, const VariantMap &) {
    return none;
  }

public:
  std::string serialize() const;

  /**
   * @brief Make object from serialized state.
   */
  static ObjectRef deserialize(const std::string &state, Context &ctx);

private:
  virtual std::string get_internal_state() const { return {}; }
  virtual void set_internal_state(std::string const &state) {}
};
} /* namespace ScriptInterface */
#endif
