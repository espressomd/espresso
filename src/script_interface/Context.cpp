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
#include "Context.hpp"
#include "ObjectState.hpp"

#include <utils/serialization/pack.hpp>

namespace ScriptInterface {
std::string Context::serialize(const ObjectHandle *o) const {
  ObjectState state;

  auto const params = o->get_parameters();
  state.params.resize(params.size());

  PackVisitor v;

  /* Pack parameters and keep track of ObjectRef parameters */
  boost::transform(params, state.params.begin(),
                   [&v](auto const &kv) -> PackedMap::value_type {
                     return {kv.first, boost::apply_visitor(v, kv.second)};
                   });

  /* Packed Object parameters */
  state.objects.resize(v.objects().size());
  boost::transform(v.objects(), state.objects.begin(), [this](auto const &kv) {
    return std::make_pair(kv.first, this->serialize(kv.second.get()));
  });

  state.name = o->name().to_string();
  state.internal_state = o->get_internal_state();

  return Utils::pack(state);
}

ObjectRef Context::deserialize(const std::string &state_, Context &ctx) {
  auto const state = Utils::unpack<ObjectState>(state_);

  std::unordered_map<ObjectId, ObjectRef> objects;
  boost::transform(state.objects, std::inserter(objects, objects.end()),
                   [&ctx](auto const &kv) {
                     return std::make_pair(kv.first,
                                           ctx.deserialize(kv.second, ctx));
                   });

  VariantMap params;
  for (auto const &kv : state.params) {
    params[kv.first] = boost::apply_visitor(UnpackVisitor(objects), kv.second);
  }

  auto o = ctx.make_shared(state.name, params);
  o->set_internal_state(state.internal_state);

  return o;
}
} // namespace ScriptInterface