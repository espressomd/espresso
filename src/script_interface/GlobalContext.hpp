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
#ifndef ESPRESSO_SCRIPT_INTERFACE_OBJECTMANAGER_HPP
#define ESPRESSO_SCRIPT_INTERFACE_OBJECTMANAGER_HPP

/** @file
 *
 *  Infrastructure to synchronize objects created on the head node
 *  with their corresponding remote copies.
 *
 *  Implementation in @ref GlobalContext.cpp.
 */

#include "Context.hpp"
#include "LocalContext.hpp"
#include "MpiCallbacks.hpp"
#include "ObjectHandle.hpp"
#include "packed_variant.hpp"

#include <utils/Factory.hpp>

#include <boost/serialization/utility.hpp>

#include <cstddef>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>

namespace ScriptInterface {

/**
 * @brief Global synchronizing context.
 *
 * Objects created in this context are synchronized
 * between multiple MPI ranks. That is, for each instance
 * created on the head node, a copy is created on all
 * other ranks. If the original copy is mutated via
 * @ref notify_set_parameter(), this change is also applied to all the
 * copies. Calls to @ref notify_call_method() are also propagated to all
 * ranks. The lifetime of the copies is tied to the original,
 * if the original copy is destroyed on the head node,
 * the remote copies are also destroyed.
 */
class GlobalContext : public Context {
  using ObjectId = std::size_t;

private:
  /* Instances on this node that are managed by the
   * head node. */
  std::unordered_map<ObjectId, ObjectRef> m_local_objects;

  std::shared_ptr<LocalContext> m_node_local_context;

  bool m_is_head_node;

private:
  Communication::CallbackHandle<ObjectId, const std::string &,
                                const PackedMap &>
      cb_make_handle;
  Communication::CallbackHandle<ObjectId, const std::string &,
                                const PackedVariant &>
      cb_set_parameter;
  Communication::CallbackHandle<ObjectId, std::string const &,
                                PackedMap const &>
      cb_call_method;
  Communication::CallbackHandle<ObjectId> cb_delete_handle;

public:
  GlobalContext(Communication::MpiCallbacks &callbacks,
                std::shared_ptr<LocalContext> node_local_context)
      : m_local_objects(), m_node_local_context(std::move(node_local_context)),
        m_is_head_node(callbacks.comm().rank() == 0),
        cb_make_handle(&callbacks,
                       [this](ObjectId id, const std::string &name,
                              const PackedMap &parameters) {
                         make_handle(id, name, parameters);
                       }),
        cb_set_parameter(&callbacks,
                         [this](ObjectId id, std::string const &name,
                                PackedVariant const &value) {
                           set_parameter(id, name, value);
                         }),
        cb_call_method(&callbacks,
                       [this](ObjectId id, std::string const &name,
                              PackedMap const &arguments) {
                         call_method(id, name, arguments);
                       }),
        cb_delete_handle(&callbacks,
                         [this](ObjectId id) { delete_handle(id); }) {}

private:
  /**
   * @brief Callback for @c cb_make_handle
   */
  void make_handle(ObjectId id, const std::string &name,
                   const PackedMap &parameters);
  /**
   * @brief Create remote instances
   *
   * @param id Internal identifier of the instance
   * @param name Class name
   * @param parameters Constructor parameters.
   */
  void remote_make_handle(ObjectId id, const std::string &name,
                          const VariantMap &parameters);

private:
  /**
   * @brief Callback for @c cb_set_parameter
   */
  void set_parameter(ObjectId id, std::string const &name,
                     PackedVariant const &value);

public:
  void notify_set_parameter(const ObjectHandle *o, std::string const &name,
                            Variant const &value) override;

private:
  /**
   * @brief Callback for @c cb_call_method
   */
  void call_method(ObjectId id, std::string const &name,
                   PackedMap const &arguments);

public:
  void notify_call_method(const ObjectHandle *o, std::string const &name,
                          VariantMap const &arguments) override;

private:
  /**
   * @brief Callback for @c cb_delete_handle
   */
  void delete_handle(ObjectId id) { m_local_objects.erase(id); }

public:
  /**
   * @brief Get a new reference counted instance of a script interface
   * object by name.
   *
   * Remote objects are automatically constructed.
   */
  std::shared_ptr<ObjectHandle>
  make_shared(std::string const &name, const VariantMap &parameters) override;

  boost::string_ref name(const ObjectHandle *o) const override;

  bool is_head_node() const override { return m_is_head_node; };
};
} // namespace ScriptInterface

#endif
