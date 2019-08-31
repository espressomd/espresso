/*
  Copyright (C) 2015-2018 The ESPResSo project

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

#ifndef SCRIPT_INTERFACE_SCRIPT_INTERFACE_BASE_HPP
#define SCRIPT_INTERFACE_SCRIPT_INTERFACE_BASE_HPP

#include <map>
#include <memory>
#include <type_traits>

#include "utils/Span.hpp"

#include <boost/utility/string_ref.hpp>

#include "Variant.hpp"

namespace ScriptInterface {
/**
 * @brief Make a Variant from argument.
 *
 * This is a convenience function, so that rather involved constructors from
 * boost::variant are not needed in the script interfaces.
 */
template <typename T> Variant make_variant(const T &x) { return Variant(x); }

/**
 * @brief Base class for generic script interfaces.
 *
 * See section @ref script_interface_howto for detailed instructions on how to
 * create derived classes.
 *
 */
class ScriptInterfaceBase : public Utils::AutoObjectId<ScriptInterfaceBase> {
public:
  enum class CreationPolicy { LOCAL, GLOBAL };

protected:
  ScriptInterfaceBase() = default;

public:
  /* Copy has unclear semantics, so it should not be allowed. */
  ScriptInterfaceBase(ScriptInterfaceBase const &) = delete;
  ScriptInterfaceBase(ScriptInterfaceBase &&) = delete;
  ScriptInterfaceBase &operator=(ScriptInterfaceBase const &) = delete;
  ScriptInterfaceBase &operator=(ScriptInterfaceBase &&) = delete;
  ~ScriptInterfaceBase() override = default;

  static std::weak_ptr<ScriptInterfaceBase> &get_instance(ObjectId id);

private:
  /* Members related to object construction, they are
     only to be used internally. */

  std::string m_name;
  CreationPolicy m_policy;

  void set_name(std::string const &name) { m_name = name; }
  void set_policy(CreationPolicy policy) { m_policy = policy; }

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
  virtual void construct(VariantMap const &params) {
    for (auto const &p : params) {
      set_parameter(p.first, p.second);
    }
  }

public:
  /**
   * @brief get current parameters.
   * @return Parameters set in class.
   */
  virtual VariantMap get_parameters() const {
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
  virtual Variant get_parameter(const std::string &name) const {
    return get_parameters().at(name);
  }

  /**
   * @brief Set single parameter.
   */
  virtual void set_parameter(const std::string &, const Variant &) {}

  /**
   * @brief Call a method on the object.
   *
   * If not overridden by the implementation, this does nothing.
   */
  virtual Variant call_method(const std::string &, const VariantMap &) {
    return none;
  }

  /**
   * @brief Get a new reference counted instance of a script interface by
   * name.
   *
   */
  static std::shared_ptr<ScriptInterfaceBase>
  make_shared(std::string const &name, CreationPolicy policy);

  /**
   * @brief Get a new reference counted instance of a script interface by
   * name, restoring the state of the object
   *
   */
  static std::shared_ptr<ScriptInterfaceBase>
  make_shared(std::string const &name, CreationPolicy policy,
              Variant const &state) {
    auto so_ptr = make_shared(name, policy);
    so_ptr->set_state(state);
    return so_ptr;
  }

  /**
   * @brief Get a new reference counted instance of a script interface by
   * type.
   *
   */
  template <typename T> std::shared_ptr<T> static make_shared() {
    std::shared_ptr<T> sp = std::make_shared<T>();

    /* id of the newly created instance */
    const auto id = sp->id();

    /* get a reference to the corresponding weak_ptr in ObjectId and update
       it with our shared ptr, so that everybody uses the same ref count.
    */
    sp->get_instance(id) = std::static_pointer_cast<ScriptInterfaceBase>(sp);

    return sp;
  }

public:
  std::string serialize() const;
  static std::shared_ptr<ScriptInterfaceBase>
  unserialize(std::string const &state);
  virtual Variant get_state() const;

protected:
  virtual void set_state(Variant const &state);
};
} /* namespace ScriptInterface */
#endif
