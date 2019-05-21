/*
 * Copyright (C) 2015-2019 The ESPResSo project
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

#include "MpiCallbacks.hpp"
#include "Variant.hpp"

#include <utils/Span.hpp>

#include <boost/optional.hpp>
#include <boost/utility/string_ref.hpp>

#include <map>
#include <memory>

namespace ScriptInterface {
/**
 * @brief Base class for generic script interfaces.
 *
 * See section @ref script_interface_howto for detailed instructions on how to
 * create derived classes.
 *
 */
class ObjectHandle {
public:
  enum class CreationPolicy { LOCAL, GLOBAL };

protected:
  ObjectHandle() = default;

public:
  /* Copy has unclear semantics, so it should not be allowed. */
  ObjectHandle(ObjectHandle const &) = delete;
  ObjectHandle &operator=(ObjectHandle const &) = delete;
  virtual ~ObjectHandle();

private:
  std::string m_name;
  CreationPolicy m_policy = CreationPolicy::LOCAL;

public:
  /**
   * @brief Name of the object.
   *
   * This is the name by which this instance was constructed.
   *
   * @return Name of the object.
   */
  std::string const &name() const { return m_name; }

  /**
   * @brief The construction policy of this instance.
   */
  CreationPolicy policy() const { return m_policy; }

  /**
   * @brief Constructor
   *
   * This function is called on object creation with user
   * provided parameters. This can be used if the SO has required parameters,
   * it represents some type that can not reasonably be default constructed,
   * or if the core implementation has to be chosen by a parameter.
   * It is guaranteed that no getter or setter functions from this interface
   * is called before construct (only name() and valid_parameters()),
   * and it is only called once.
   *
   * The default implementation just calls set_parameter for every parameter.
   *
   * @param params The parameters to the constructor. Only parameters that
   *               are valid for a default-constructed object are valid.
   */
  void construct(VariantMap const &params, CreationPolicy policy,
                 const std::string &name);

private:
  virtual void do_construct(VariantMap const &params) {
    for (auto const &p : params) {
      do_set_parameter(p.first, p.second);
    }
  }

public:
  /**
   * @brief get current parameters.
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
   * @brief Get required and optional parameters for class
   *
   * Get required and optional parameters for class.
   *
   * @return Expected parameters.
   */
  virtual Utils::Span<const boost::string_ref> valid_parameters() const {
    return {};
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
  virtual void do_set_parameter(const std::string &, const Variant &) {}

public:
  /**
   * @brief Call a method on the object.
   *
   * If not overridden by the implementation, this does nothing.
   */
  Variant call_method(const std::string &name, const VariantMap &params);

private:
  virtual Variant do_call_method(const std::string &, const VariantMap &) {
    return none;
  }

public:
  static void initialize(Communication::MpiCallbacks &cb);

  /**
   * @brief Get a new reference counted instance of a script interface by
   * name.
   *
   */
  static std::shared_ptr<ObjectHandle>
  make_shared(std::string const &name, CreationPolicy policy,
              const VariantMap &parameters = {});

public:
  std::string serialize() const;
  static std::shared_ptr<ObjectHandle> unserialize(std::string const &state);
};
} /* namespace ScriptInterface */
#endif
