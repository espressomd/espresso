/*
 * Copyright (C) 2020 The ESPResSo project
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
#ifndef ESPRESSO_CONTEXT_HPP
#define ESPRESSO_CONTEXT_HPP

#include "ObjectHandle.hpp"
#include "Variant.hpp"

#include <boost/utility/string_ref.hpp>

#include <memory>

namespace ScriptInterface {
class Context : public std::enable_shared_from_this<Context> {
public:
  /**
   * @brief Call method on remote instances
   *
   * @param self Internal identified of the instance
   * @param name Name of the method to call
   * @param arguments Arguments to the call
   */
  virtual void notify_call_method(const ObjectHandle *self,
                                  std::string const &name,
                                  VariantMap const &arguments) = 0;

  /**
   * @brief Set a parameter on remote instances
   *
   * @param self Internal identifier of the instance to be modified
   * @param name Name of the parameter to change
   * @param value Value to set it to
   */
  virtual void notify_set_parameter(const ObjectHandle *self,
                                    std::string const &name,
                                    Variant const &value) = 0;

  /**
   * @brief Get a new reference counted instance of a script interface by
   * name.
   *
   */
  virtual std::shared_ptr<ObjectHandle>
  make_shared(std::string const &name, const VariantMap &parameters) = 0;

  /**
   * @brief String representation of the state of an object.
   */
  std::string serialize(const ObjectHandle *o) const;

  /**
   * @brief Make object from serialized state.
   */
  ObjectRef deserialize(std::string const &state_);

protected:
  void set_manager(ObjectHandle *o) { o->m_manager = this->shared_from_this(); }
  void set_name(ObjectHandle *o, boost::string_ref name) const {
    o->m_name = name;
  }

public:
  boost::string_ref name(const ObjectHandle *o) const { return o->name(); }

public:
  virtual ~Context() = default;
};
} // namespace ScriptInterface
#endif // ESPRESSO_CONTEXT_HPP
