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
#include "ObjectManager.hpp"

#include "GlobalContext.hpp"
#include "LocalContext.hpp"

#include <utils/serialization/pack.hpp>

namespace ScriptInterface {
std::shared_ptr<ObjectHandle>
ObjectManager::make_shared(CreationPolicy policy, std::string const &name,
                           const VariantMap &parameters) {
  return context(policy)->make_shared(name, parameters);
}

std::shared_ptr<ObjectHandle>
ObjectManager::deserialize(std::string const &state_) {
  auto const state =
      Utils::unpack<std::pair<CreationPolicy, std::string>>(state_);

  return context(state.first)->deserialize(state.second);
}

std::string ObjectManager::serialize(const ObjectHandle *o) const {
  /* We treat objects without a context as local. */
  auto ctx = o->context() ? o->context() : m_local_context.get();

  return Utils::pack(std::make_pair(policy(ctx), ctx->serialize(o)));
}

ObjectManager::ObjectManager(Communication::MpiCallbacks &callbacks,
                             const Utils::Factory<ObjectHandle> &factory) {
  auto local_context = std::make_shared<LocalContext>(factory);

  /* If there is only one node, we can treat all objects as local, and thus
   * never invoke any callback. */
  m_global_context =
      (callbacks.comm().size() > 1)
          ? std::make_shared<GlobalContext>(callbacks, local_context)
          : std::static_pointer_cast<Context>(local_context);

  m_local_context = std::move(local_context);
}
} // namespace ScriptInterface
