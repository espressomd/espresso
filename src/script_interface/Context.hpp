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
#ifndef ESPRESSO_SCRIPT_INTERFACE_CONTEXT_HPP
#define ESPRESSO_SCRIPT_INTERFACE_CONTEXT_HPP

/** @file
 *
 *  @ref ScriptInterface::Context decorates @ref ScriptInterface::ObjectHandle
 *  objects with a context: a creation policy (local object, local object
 *  with remote copies) and a communication facility to synchronize an object
 *  on the head node with remote copies (serialization, callback mechanism).
 */

#include "ObjectHandle.hpp"
#include "Variant.hpp"

#include <boost/utility/string_ref.hpp>

#include <memory>
#include <string>

namespace ScriptInterface {
/**
 * @brief Context of an object handle.
 *
 * Each instance of @ref ObjectHandle can have an
 * attached context, which can e.g. synchronize
 * distributed copies of the instance. The context
 * does also provide facilities for serializing
 * Objects into a string representation.
 */
class Context : public std::enable_shared_from_this<Context> {
public:
  /**
   * @brief Call method on remote instances
   *
   * @param self Internal identifier of the instance
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
   * Objects created thru a Context get shared ownership of that context,
   * e.g. the lifetime of the context is at least as long as the objects
   * created by it. Therefore the object can always assume that the context
   * is present.
   */
  virtual std::shared_ptr<ObjectHandle>
  make_shared(std::string const &name, const VariantMap &parameters) = 0;

protected:
  /**
   * @brief Set the context of an object to this.
   *
   * @param o Object to set the context for.
   */
  void set_context(ObjectHandle *o) { o->m_context = this->shared_from_this(); }

public:
  /**
   *  @brief Get the class name for an ObjectHandle instance.
   *
   *  This returns the name by which the object can be created.
   */
  virtual boost::string_ref name(const ObjectHandle *o) const = 0;

  virtual bool is_head_node() const = 0;

  virtual ~Context() = default;
};
} // namespace ScriptInterface
#endif // ESPRESSO_CONTEXT_HPP
