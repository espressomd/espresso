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
#include "ScriptInterface.hpp"
#include "Serializer.hpp"
#include "pack.hpp"
#include "PackedVariant.hpp"

#include <utils/serialization/pack.hpp>

#include <boost/range/algorithm/for_each.hpp>
#include <boost/serialization/utility.hpp>

namespace ScriptInterface {
namespace {
Communication::MpiCallbacks *m_callbacks = nullptr;

void make_remote_handle(ObjectId id, const std::string &name,
                        const PackedMap &parameters) {
  local_objects[id] = ObjectHandle::make_shared(
      name, ObjectHandle::CreationPolicy::LOCAL, unpack(parameters));
}

void remote_set_parameter(ObjectId id, std::string const &name,
                          PackedVariant const &value) {
  local_objects.at(id)->set_parameter(name, unpack(value));
}

void remote_call_method(ObjectId id, std::string const &name,
                        PackedMap const &arguments) {
  local_objects.at(id)->call_method(name, unpack(arguments));
}

void delete_remote_handle(ObjectId id) { local_objects.erase(id); }

REGISTER_CALLBACK(make_remote_handle)
REGISTER_CALLBACK(remote_set_parameter)
REGISTER_CALLBACK(remote_call_method)
REGISTER_CALLBACK(delete_remote_handle)
} // namespace

Utils::Factory<ObjectHandle> factory;

std::shared_ptr<ObjectHandle>
ObjectHandle::make_shared(std::string const &name, CreationPolicy policy,
                          const VariantMap &parameters) {
  std::shared_ptr<ObjectHandle> sp = factory.make(name);

  sp->construct(parameters, policy, name);

  return sp;
}

/**
 * @brief Returns a binary representation of the state often
 *        the instance, as returned by get_state().
 */
std::string ObjectHandle::serialize() const {
  //return Utils::pack(Serializer{}(this));
  return {};
}

/**
 * @brief Creates a new instance from a binary state,
 *        as returned by serialize().
 */
std::shared_ptr<ObjectHandle>
ObjectHandle::unserialize(std::string const &state) {
  return {};
}

void ObjectHandle::construct(VariantMap const &params, CreationPolicy policy,
                             const std::string &name) {
  m_name = name;
  m_policy = policy;

  if (m_policy == CreationPolicy::GLOBAL) {
    m_callbacks->call(make_remote_handle, object_id(this), name, pack(params));
  }

  this->do_construct(params);
}

void ObjectHandle::set_parameter(const std::string &name,
                                 const Variant &value) {
  if (m_policy == CreationPolicy::GLOBAL) {
    m_callbacks->call(remote_set_parameter, object_id(this), name, pack(value));
  }

  this->do_set_parameter(name, value);
}

Variant ObjectHandle::call_method(const std::string &name,
                                  const VariantMap &params) {
  if (m_policy == CreationPolicy::GLOBAL) {
    m_callbacks->call(remote_call_method, object_id(this), name, pack(params));
  }

  return this->do_call_method(name, params);
}

ObjectHandle::~ObjectHandle() {
  if (m_policy == CreationPolicy::GLOBAL) {
    m_callbacks->call(delete_remote_handle, object_id(this));
  }
}

void ObjectHandle::initialize(::Communication::MpiCallbacks &cb) {
  m_callbacks = &cb;
}

Variant ObjectHandle::get_state() const {
  std::vector<Variant> state;

  auto params = this->get_parameters();
  state.reserve(3 + params.size());

  state.push_back(static_cast<int>(m_policy));
  state.push_back(m_name);

  for (auto const &p : params) {
    state.push_back(std::vector<Variant>{
        {p.first, boost::apply_visitor(Serializer{}, p.second)}});
  }

  return state;
}

void ObjectHandle::set_state(Variant const &state) {
  using boost::make_iterator_range;
  using boost::get;
  using std::vector;

  auto const& state_ = get<vector<Variant>>(state);
  auto const policy = CreationPolicy(get<int>(state_.at(0)));
  auto const& name = get<std::string>(state_.at(1));

  UnSerializer u;
  VariantMap params;

  for (auto const &v : make_iterator_range(state_.begin() + 2, state_.end())) {
    auto const &p = get<vector<Variant>>(v);
    params[get<std::string>(p.at(0))] = boost::apply_visitor(u, p.at(1));
  }

  this->construct(params, policy, name);
}
} /* namespace ScriptInterface */
