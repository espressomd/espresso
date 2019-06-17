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
#include "PackedVariant.hpp"
#include "ScriptInterface.hpp"

#include <utils/serialization/pack.hpp>

#include <boost/range/algorithm/for_each.hpp>
#include <boost/serialization/utility.hpp>

namespace ScriptInterface {
namespace {
Communication::MpiCallbacks *m_callbacks = nullptr;
std::unordered_map<ObjectId, ObjectRef> local_objects;

/**
 * @brief Callback for creating a local instance
 *
 * @param id Internal identifier of the instance
 * @param name Class name
 * @param parameters Constructor parameters.
 */
void make_remote_handle(ObjectId id, const std::string &name,
                        const PackedMap &parameters) {
  local_objects[id] =
      ObjectHandle::make_shared(name, ObjectHandle::CreationPolicy::LOCAL,
                                unpack(parameters, local_objects));
}

/**
 * @brief Callback for setting a parameter on an instance
 *
 * @param id Internal identifier of the instance to be modified
 * @param name Name of the parameter to change
 * @param value Value to set it to
 */
void remote_set_parameter(ObjectId id, std::string const &name,
                          PackedVariant const &value) {
  local_objects.at(id)->set_parameter(name, unpack(value, local_objects));
}

/**
 * @brief Callback for calling a method on an instance
 *
 * @param id Internal identified of the instance
 * @param name Name of the method to call
 * @param arguments Arguments to the call
 */
void remote_call_method(ObjectId id, std::string const &name,
                        PackedMap const &arguments) {
  local_objects.at(id)->call_method(name, unpack(arguments, local_objects));
}

/**
 * @brief Callback for deleting an instance
 *
 * @param id Internal identified of the instance
 */
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
  auto sp = factory.make(name);

  sp->construct(parameters, policy, name);

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

  template <class Archive> void serialize(Archive &ar, long int) {
    ar &name &policy &params &objects;
  }
};

/**
 * @brief Returns a binary representation of the state often
 *        the instance, as returned by get_state().
 */
std::string ObjectHandle::serialize() const {
  ObjectState state{
    name(), policy(), {}, {}
  };

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

/**
 * @brief Creates a new instance from a binary state,
 *        as returned by serialize().
 */
std::shared_ptr<ObjectHandle>
ObjectHandle::unserialize(std::string const &state_) {
  auto state = Utils::unpack<ObjectState>(state_);

  std::unordered_map<ObjectId, ObjectRef> objects;
  boost::transform(state.objects, std::inserter(objects, objects.end()), [](auto const& kv) {
    return std::make_pair(kv.first, unserialize(kv.second));
  });

  VariantMap params;
  for(auto const&kv: state.params) {
    params[kv.first] = boost::apply_visitor(UnpackVisitor(objects), kv.second);
  }

  return make_shared(state.name, state.policy, params);
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

PackedVariant ObjectHandle::get_state() const {}
void ObjectHandle::set_state(Variant const &state) {}
} /* namespace ScriptInterface */
