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
#include "Variant.hpp"

#include <utils/Span.hpp>

#include <boost/mpi/communicator.hpp>
#include <boost/utility/string_ref.hpp>

namespace ScriptInterface {
class ObjectManager;

/**
 * @brief Base class for generic script interfaces.
 *
 * See section @ref script_interface_howto for detailed instructions on how to
 * create derived classes.
 *
 */
class ObjectHandle {
protected:
  ObjectHandle() = default;

public:
  /* Copy has unclear semantics, so it should not be allowed. */
  ObjectHandle(ObjectHandle const &) = delete;
  ObjectHandle &operator=(ObjectHandle const &) = delete;
  virtual ~ObjectHandle();

private:
  friend class ObjectManager;
  std::shared_ptr<ObjectManager> m_manager = {};
  std::string m_name;
  CreationPolicy m_policy = CreationPolicy::LOCAL;

public:
  /**
   * @brief Responsible manager.
   */
  ObjectManager *manager() const {
    Utils::print(__PRETTY_FUNCTION__, boost::mpi::communicator().rank());
    auto const ptr = m_manager.get();
    return ptr; }

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

private:
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
  /**
   * @brief Local implementation of @function set_parameter
   */
  virtual void do_set_parameter(const std::string &, const Variant &) {}

public:
  /**
   * @brief Call a method on the object.
   */
  Variant call_method(const std::string &name, const VariantMap &params);

private:
  /**
   * @brief Local implementation of @do_call_method
   *
   * If not overridden by the implementation, this does nothing.
   */
  virtual Variant do_call_method(const std::string &, const VariantMap &) {
    return none;
  }

private:
  virtual std::string get_internal_state() const { return {}; }
  virtual void set_internal_state(std::string const &state) {}

  /**
   * @brief Call from the destructor of the Handle.
   *
   * This can use to customize the deletion of objects, e.g.
   * to change when the remote objects are deleted.
   * The default implementation just deletes all remote instances.
   *
   */
  virtual void do_destroy() { delete_remote(); }

protected:
  /**
   * @brief Delete remote instance of this object.
   *
   * Explicitly delete instances on other nodes, if any.
   * This can only be called once.
   */
  void delete_remote();
};
} /* namespace ScriptInterface */
#endif
