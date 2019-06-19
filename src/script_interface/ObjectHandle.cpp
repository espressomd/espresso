/*
  Copyright (C) 2010-2018 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
  Max-Planck-Institute for Polymer Research, Theory Group

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

#include "ObjectHandle.hpp"

#include "MpiCallbacks.hpp"
#include "ObjectManager.hpp"
#include "ScriptInterface.hpp"
#include "PackedVariant.hpp"

#include <utils/serialization/pack.hpp>

#include <boost/serialization/utility.hpp>
#include <boost/range/algorithm/for_each.hpp>
#include <boost/range/algorithm/transform.hpp>

namespace ScriptInterface {
namespace {
std::unique_ptr<ObjectManager> m_om;
}

Utils::Factory<ObjectHandle> factory;

std::shared_ptr<ObjectHandle>
ObjectHandle::make_shared(std::string const &name, CreationPolicy policy,
                          const VariantMap &parameters) {
  auto sp = factory.make(name);

  sp->m_name = name;
  sp->m_policy = policy;

  if (sp->m_policy == CreationPolicy::GLOBAL) {
    m_om->remote_make_handle(object_id(sp.get()), name, parameters);
  }

  sp->do_construct(parameters);

  return sp;
}

/**
 * @brief State of an object ready for serialization.
 *
 * This specifies the internal serialization format and
 * should not be used outside of the class.
 */
struct ObjectState {
  std::string name;
  ObjectHandle::CreationPolicy policy;
  PackedMap params;
  std::vector<std::pair<ObjectId, std::string>> objects;
  std::string internal_state;

  template <class Archive> void serialize(Archive &ar, long int) {
    ar &name &policy &params &objects &internal_state;
  }
};

std::string ObjectHandle::serialize() const {
  ObjectState state{name(), policy(), {}, {}, get_internal_state()};

  auto const params = get_parameters();
  state.params.resize(params.size());

  PackVisitor v;

  /* Pack parameters and keep track of ObjectRef parameters */
  boost::transform(params, state.params.begin(),
                   [&v](auto const &kv) -> PackedMap::value_type {
                     return {kv.first, boost::apply_visitor(v, kv.second)};
                   });

  /* Packed Object parameters */
  state.objects.resize(v.objects().size());
  boost::transform(v.objects(), state.objects.begin(), [](auto const &kv) {
    return std::make_pair(kv.first, kv.second->serialize());
  });

  return Utils::pack(state);
}

std::shared_ptr<ObjectHandle>
ObjectHandle::unserialize(std::string const &state_) {
  auto state = Utils::unpack<ObjectState>(state_);

  std::unordered_map<ObjectId, ObjectRef> objects;
  boost::transform(state.objects, std::inserter(objects, objects.end()),
                   [](auto const &kv) {
                     return std::make_pair(kv.first, unserialize(kv.second));
                   });

  VariantMap params;
  for (auto const &kv : state.params) {
    params[kv.first] = boost::apply_visitor(UnpackVisitor(objects), kv.second);
  }

  auto so = make_shared(state.name, state.policy, params);
  so->set_internal_state(state.internal_state);

  return so;
}

void ObjectHandle::set_parameter(const std::string &name,
                                 const Variant &value) {
  if (m_policy == CreationPolicy::GLOBAL) {
    assert(m_om), m_om->remote_set_parameter(object_id(this), name, value);
  }

  this->do_set_parameter(name, value);
}

Variant ObjectHandle::call_method(const std::string &name,
                                  const VariantMap &params) {
  if (m_policy == CreationPolicy::GLOBAL) {
    m_om->remote_call_method(object_id(this), name, params);
  }

  return this->do_call_method(name, params);
}

void ObjectHandle::delete_remote() {
  if (m_policy == CreationPolicy::GLOBAL) {
    m_om->remote_delete_handle(object_id(this));
  }
}

ObjectHandle::~ObjectHandle() { this->do_destroy(); }

void ObjectHandle::initialize(::Communication::MpiCallbacks &cb) {
  m_om = std::make_unique<ObjectManager>(&cb);
}
} /* namespace ScriptInterface */
