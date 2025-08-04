/*
 * Copyright (C) 2010-2022 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 * Max-Planck-Institute for Polymer Research, Theory Group
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

#include "ObjectHandle.hpp"
#include "Context.hpp"
#include "ObjectState.hpp"
#include "packed_variant.hpp"

#include <config/config.hpp>

#include <instrumentation/fe_trap.hpp>

#include <utils/serialization/pack.hpp>

#include <algorithm>
#include <iterator>
#include <memory>
#include <string>
#include <string_view>
#include <unordered_map>
#include <utility>

namespace ScriptInterface {
void ObjectHandle::set_parameter(const std::string &name,
                                 const Variant &value) {
  if (m_context)
    m_context->notify_set_parameter(this, name, value);

#ifdef FPE
  auto const trap = fe_trap::make_shared_scoped();
#endif
  this->do_set_parameter(name, value);
}

Variant ObjectHandle::call_method(const std::string &name,
                                  const VariantMap &params) {
  if (m_context)
    m_context->notify_call_method(this, name, params);

#ifdef FPE
  auto const trap = fe_trap::make_shared_scoped();
#endif
  return this->do_call_method(name, params);
}

std::string ObjectHandle::serialize() const {
  ObjectState state;

  auto const params = serialize_parameters();
  state.params.reserve(params.size());

  PackVisitor visit;

  /* Pack parameters and keep track of ObjectRef parameters */
  std::ranges::transform(params, std::back_inserter(state.params),
                         [&visit](auto const &kv) -> PackedMap::value_type {
                           auto const &[name, value] = kv;
                           return {name, boost::apply_visitor(visit, value)};
                         });

  /* Packed Object parameters */
  state.objects.reserve(visit.objects().size());
  std::ranges::transform(visit.objects(), std::back_inserter(state.objects),
                         [](auto const &kv) {
                           auto const &[name, obj] = kv;
                           return std::make_pair(name, obj->serialize());
                         });

  state.name = name();
  state.internal_state = get_internal_state();

  return Utils::pack(state);
}

ObjectRef ObjectHandle::deserialize(const std::string &packed_state,
                                    Context &ctx) {
  auto const state = Utils::unpack<ObjectState>(packed_state);

  std::unordered_map<ObjectId, ObjectRef> objects;
  std::ranges::transform(state.objects, std::inserter(objects, objects.end()),
                         [&ctx](auto const &kv) {
                           auto const &[name, buf] = kv;
                           return std::make_pair(name, deserialize(buf, ctx));
                         });

  VariantMap params;
  for (auto const &[name, variant] : state.params) {
    params[name] = boost::apply_visitor(UnpackVisitor(objects), variant);
  }

  auto o = ctx.make_shared(state.name, params);
  o->set_internal_state(state.internal_state);

  return o;
}

std::string_view ObjectHandle::name() const {
  return context() ? context()->name(this) : std::string_view{};
}

} /* namespace ScriptInterface */
