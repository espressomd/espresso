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
#ifndef ESPRESSO_CONTEXTMANAGER_HPP
#define ESPRESSO_CONTEXTMANAGER_HPP

/** @file
 *
 *  @ref ScriptInterface::ContextManager manages object creation with policies
 *  @ref ScriptInterface::ContextManager::CreationPolicy "CreationPolicy".
 *  Object creation is delegated to @ref ScriptInterface::GlobalContext and
 *  @ref ScriptInterface::LocalContext. @ref ScriptInterface::ContextManager
 *  serves as their public interface. If there is only 1 MPI rank, no
 *  communication takes place and all objects are created locally via
 *  @ref ScriptInterface::LocalContext, including those with policy
 *  @ref ScriptInterface::ContextManager::CreationPolicy::GLOBAL "GLOBAL".
 *
 *  Implementation in @ref ContextManager.cpp.
 */

#include "Context.hpp"
#include "Variant.hpp"

#include "core/MpiCallbacks.hpp"

#include <utils/Factory.hpp>

#include <cassert>
#include <memory>
#include <stdexcept>
#include <string>

namespace ScriptInterface {

/**
 * @brief Manage object contexts.
 *
 * This owns object contexts and allows for
 * creation and serialization of objects preserving
 * their context.
 */
class ContextManager {
  std::shared_ptr<Context> m_local_context;
  std::shared_ptr<Context> m_global_context;

public:
  /** Labels for context */
  enum class CreationPolicy {
    /** Corresponding to @c LocalContext */
    LOCAL,
    /** Corresponding to @c GlobalContext */
    GLOBAL
  };

  ContextManager(Communication::MpiCallbacks &callbacks,
                 const Utils::Factory<ObjectHandle> &factory);

  /**
   * @brief Get a new reference counted instance of a script interface by
   * name.
   */
  std::shared_ptr<ObjectHandle> make_shared(CreationPolicy policy,
                                            std::string const &name,
                                            const VariantMap &parameters);

  /**
   * @brief Get a new reference counted instance of a script interface from
   * a serialized state.
   */
  std::shared_ptr<ObjectHandle> deserialize(std::string const &state_);

  /**
   * @brief Serialize a script interface object into a binary representation.
   */
  std::string serialize(const ObjectHandle *o) const;

private:
  /**
   * @brief Map policy to context.
   *
   * Inverse of policy.
   */
  Context *context(CreationPolicy policy) const {
    switch (policy) {
    case CreationPolicy::LOCAL:
      return assert(m_local_context), m_local_context.get();
    case CreationPolicy::GLOBAL:
      return assert(m_global_context), m_global_context.get();
    default:
      throw std::runtime_error("Unknown context type.");
    }
  }

  /**
   * @brief Map context to policy.
   *
   * Inverse of context.
   */
  CreationPolicy policy(Context *c) const {
    if (c == m_local_context.get()) {
      return CreationPolicy::LOCAL;
    }
    if (c == m_global_context.get()) {
      return CreationPolicy::GLOBAL;
    }

    throw std::runtime_error("Invalid context.");
  }
};
} // namespace ScriptInterface

#endif // ESPRESSO_CONTEXTMANAGER_HPP
