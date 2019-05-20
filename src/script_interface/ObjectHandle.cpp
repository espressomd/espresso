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

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/stream.hpp>
#include <boost/serialization/utility.hpp>

#include <sstream>

namespace ScriptInterface {
namespace {
std::unordered_map<ObjectId, std::shared_ptr<ObjectHandle>> local_objects;

using PackedVariant = boost::make_recursive_variant<
    None, bool, int, size_t, double, std::string, std::vector<int>, std::vector<double>,
    ObjectId, std::vector<boost::recursive_variant_>, Utils::Vector2d,
    Utils::Vector3d, Utils::Vector4d>::type;

using PackedMap = std::vector<std::pair<std::string, PackedVariant>>;

template <class D, class V, class R>
struct recursive_visitor : boost::static_visitor<R> {
  template <class T> R operator()(T &&val) const {
    return std::forward<T>(val);
  }

  R operator()(const std::vector<V> &vec) const {
    std::vector<R> ret(vec.size());

    boost::transform(vec, ret.begin(),
                     [visitor = static_cast<const D *>(this)](const V &v) {
                       return boost::apply_visitor(*visitor, v);
                     });

    return ret;
  }
};

struct VariantToTransport
    : recursive_visitor<VariantToTransport, Variant, PackedVariant> {
  using recursive_visitor<VariantToTransport, Variant, PackedVariant>::
  operator();

  PackedVariant operator()(const ObjectRef &so_ptr) const {
    return so_ptr->id();
  }
};

struct TransportToVariant
    : recursive_visitor<TransportToVariant, PackedVariant, Variant> {
  using recursive_visitor<TransportToVariant, PackedVariant, Variant>::
  operator();

  Variant operator()(const ObjectId &id) const { return local_objects.at(id); }
};

PackedVariant pack(const Variant &v) {
  return boost::apply_visitor(VariantToTransport{}, v);
}

Variant unpack(const PackedVariant &v) {
  return boost::apply_visitor(TransportToVariant{}, v);
}

PackedMap pack(const VariantMap &v) {
  std::vector<std::pair<std::string, PackedVariant>> ret(v.size());

  boost::transform(v, ret.begin(), [](auto const &kv) {
    return std::pair<std::string, PackedVariant>{kv.first, pack(kv.second)};
  });

  return ret;
}

VariantMap unpack(const PackedMap &v) {
  VariantMap ret;

  boost::transform(v, std::inserter(ret, ret.end()), [](auto const &kv) {
    return std::pair<std::string, Variant>{kv.first, unpack(kv.second)};
  });

  return ret;
}

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

void delete_remote_handle(ObjectId id) {
  local_objects.erase(id);
}

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
std::string ObjectHandle::serialize() const { return {}; }

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
    m_callbacks->call(make_remote_handle, id(), name, pack(params));
  }

  this->do_construct(params);
}

void ObjectHandle::set_parameter(const std::string &name,
                                 const Variant &value) {
  if (m_policy == CreationPolicy::GLOBAL) {
    m_callbacks->call(remote_set_parameter, id(), name, pack(value));
  }

  this->do_set_parameter(name, value);
}

Variant ObjectHandle::call_method(const std::string &name,
                                  const VariantMap &params) {
  if (m_policy == CreationPolicy::GLOBAL) {
    m_callbacks->call(remote_call_method, id(), name, pack(params));
  }

  return this->do_call_method(name, params);
}

ObjectHandle::~ObjectHandle() {
  if (m_policy == CreationPolicy::GLOBAL) {
    m_callbacks->call(delete_remote_handle, id());
  }
}

void ObjectHandle::initialize(::Communication::MpiCallbacks &cb) {
  m_callbacks = &cb;
}
} /* namespace ScriptInterface */
