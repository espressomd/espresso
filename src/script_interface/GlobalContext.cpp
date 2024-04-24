/*
 * Copyright (C) 2020-2022 The ESPResSo project
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

/** @file
 *
 *  Infrastructure to synchronize objects created on the head node with
 *  their corresponding remote copies. Methods of script interface
 *  classes may throw exceptions of type @ref ScriptInterface::Exception.
 *  These exceptions will halt the flow of the program on the head node.
 *  The same exceptions will be thrown in the remote copies but will
 *  be silenced, since they are redundant. Other types of exceptions
 *  are not silenced.
 *
 *  Implementation of @ref GlobalContext.hpp.
 */

#include "GlobalContext.hpp"
#include "Exception.hpp"
#include "ObjectHandle.hpp"
#include "packed_variant.hpp"

#include <utils/serialization/pack.hpp>

#include <cassert>
#include <stdexcept>
#include <string>
#include <utility>

namespace ScriptInterface {
void GlobalContext::make_handle(ObjectId id, const std::string &name,
                                const PackedMap &parameters) {
  try {
    ObjectRef so = m_node_local_context->make_shared(
        name, unpack(parameters, m_local_objects));

    m_local_objects[id] = std::move(so);
  } catch (Exception const &) { // NOLINT(bugprone-empty-catch)
  }
}

void GlobalContext::remote_make_handle(ObjectId id, const std::string &name,
                                       const VariantMap &parameters) {
  cb_make_handle(id, name, pack(parameters));
}

void GlobalContext::set_parameter(ObjectId id, std::string const &name,
                                  PackedVariant const &value) {
  try {
    m_local_objects.at(id)->set_parameter(name, unpack(value, m_local_objects));
  } catch (Exception const &) { // NOLINT(bugprone-empty-catch)
  }
}

void GlobalContext::notify_set_parameter(const ObjectHandle *o,
                                         std::string const &name,
                                         Variant const &value) {
  cb_set_parameter(object_id(o), name, pack(value));
}

void GlobalContext::call_method(ObjectId id, std::string const &name,
                                PackedMap const &arguments) {
  try {
    m_local_objects.at(id)->call_method(name,
                                        unpack(arguments, m_local_objects));
  } catch (Exception const &) { // NOLINT(bugprone-empty-catch)
  }
}

void GlobalContext::notify_call_method(const ObjectHandle *o,
                                       std::string const &name,
                                       VariantMap const &arguments) {
  cb_call_method(object_id(o), name, pack(arguments));
}

std::shared_ptr<ObjectHandle>
GlobalContext::make_shared_local(std::string const &name,
                                 VariantMap const &parameters) {
  auto sp = m_node_local_context->factory().make(name);
  set_context(sp.get());

  sp->construct(parameters);

  return sp;
}

std::shared_ptr<ObjectHandle>
GlobalContext::make_shared(std::string const &name,
                           const VariantMap &parameters) {
  assert(is_head_node());
  auto sp = m_node_local_context->factory().make(name);
  set_context(sp.get());

  auto const id = object_id(sp.get());
  remote_make_handle(id, name, parameters);

  sp->construct(parameters);

  return {sp.release(),
          /* Custom deleter, we keep the corresponding global context,
           * as well as the original deleter for the object. */
          [global_context = this, deleter = sp.get_deleter()](ObjectHandle *o) {
            /* Tell the other nodes before invoking the destructor, this is
             * required
             * to have synchronous destructors, which is needed by some client
             * code. */
            global_context->cb_delete_handle(object_id(o));

            /* Locally destroy the object. */
            deleter(o);
          }};
}

boost::string_ref GlobalContext::name(const ObjectHandle *o) const {
  assert(o);

  return m_node_local_context->factory().type_name(*o);
}
} // namespace ScriptInterface