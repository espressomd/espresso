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
#include <boost/serialization/variant.hpp>

#include <sstream>

namespace ScriptInterface {
std::unordered_map<ObjectId, std::weak_ptr<ObjectHandle>> local_objects;

namespace detail {
using TransportVariant = boost::make_recursive_variant<
    None, bool, int, size_t, double, std::string, std::vector<int>,
    std::vector<double>, ObjectId, std::vector<boost::recursive_variant_>,
    Utils::Vector2d, Utils::Vector3d, Utils::Vector4d>::type;

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
    : recursive_visitor<VariantToTransport, Variant, TransportVariant> {
  using recursive_visitor<VariantToTransport, Variant, TransportVariant>::
  operator();

  TransportVariant operator()(const ObjectRef &so_ptr) const {
    return so_ptr->id();
  }
};

struct TransportToVariant
    : recursive_visitor<TransportToVariant, TransportVariant, Variant> {
  using recursive_visitor<TransportToVariant, TransportVariant, Variant>::
  operator();

  Variant operator()(const ObjectId &id) const {
    return local_objects.at(id).lock();
  }
};

TransportVariant pack(const Variant &v) {
  return boost::apply_visitor(VariantToTransport{}, v);
}

Variant unpack(const TransportVariant &v) {
  return boost::apply_visitor(TransportToVariant{}, v);
}

std::vector<std::pair<std::string, TransportVariant>>
pack(const VariantMap &v) {
  std::vector<std::pair<std::string, TransportVariant>> ret(v.size());

  boost::transform(v, ret.begin(), [](auto const &kv) {
    return std::pair<std::string, TransportVariant>{kv.first, pack(kv.second)};
  });

  return ret;
}

VariantMap
unpack(const std::vector<std::pair<std::string, TransportVariant>> &v) {
  VariantMap ret;

  boost::transform(v, std::inserter(ret, ret.end()), [](auto const &kv) {
    return std::pair<std::string, Variant>{kv.first, unpack(kv.second)};
  });

  return ret;
}

struct CallbackAction {
  struct Construct {
    ObjectId id;
    std::string name;
    std::vector<std::pair<std::string, TransportVariant>> parameters;

    template <class Archive> void serialize(Archive &ar, long int) {
      ar &id &name &parameters;
    }
  };
  struct SetParameter {
    std::string name;
    TransportVariant value;

    template <class Archive> void serialize(Archive &ar, long int) {
      ar &name &value;
    }
  };
  struct CallMethod {
    std::string name;
    std::vector<std::pair<std::string, TransportVariant>> arguments;

    template <class Archive> void serialize(Archive &ar, long int) {
      ar &name &arguments;
    }
  };

  boost::variant<Construct, SetParameter, CallMethod> value;

  template <class Archive> void serialize(Archive &ar, long int) { ar &value; }
};
} // namespace detail

namespace {
Communication::MpiCallbacks *m_callbacks = nullptr;

class RemoteObjectHandle {
private:
  detail::Callback m_callback_id;

  const auto &comm() const { return m_callback_id.cb()->comm(); }
  std::shared_ptr<ObjectHandle> m_p;

public:
  RemoteObjectHandle(Communication::MpiCallbacks *cb)
      : m_callback_id(cb, [this](detail::CallbackAction a) { mpi_slave(a); }) {}

  struct CallbackVisitor {
    using CallbackAction = detail::CallbackAction;

    std::shared_ptr<ObjectHandle> &o;

    void operator()(const CallbackAction::Construct &ctor) const {
      o = ObjectHandle::make_shared(ctor.name,
                                    ObjectHandle::CreationPolicy::LOCAL,
                                    detail::unpack(ctor.parameters));
      local_objects[ctor.id] = o;
    }
    void operator()(const CallbackAction::SetParameter &param) const {
      assert(o), o->set_parameter(param.name, detail::unpack(param.value));
    }
    void operator()(const CallbackAction::CallMethod &method) const {
      assert(o);
      (void)o->call_method(method.name, detail::unpack(method.arguments));
    }
  };

  void mpi_slave(detail::CallbackAction a) {
    boost::apply_visitor(CallbackVisitor{this->m_p}, a.value);
  }
};
} // namespace
} // namespace ScriptInterface

static void make_remote_handle() {
  using namespace ScriptInterface;
  new RemoteObjectHandle(m_callbacks);
}

static void delete_remote_handle() {
  using namespace ScriptInterface;
  /* TODO: Implement */
}

REGISTER_CALLBACK(make_remote_handle)
REGISTER_CALLBACK(delete_remote_handle)

namespace ScriptInterface {
Utils::Factory<ObjectHandle> factory;

std::shared_ptr<ObjectHandle>
ObjectHandle::make_shared(std::string const &name, CreationPolicy policy,
                          const VariantMap &parameters) {
  std::shared_ptr<ObjectHandle> sp = factory.make(name);

  sp->construct(parameters, policy, name);

  return sp;
}

std::weak_ptr<ObjectHandle> &ObjectHandle::get_instance(ObjectId id) {
  return Utils::AutoObjectId<ObjectHandle>::get_instance(id);
}

/* Checkpointing functions. */

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
  using detail::CallbackAction;
  using Construct = CallbackAction::Construct;

  m_name = name;
  m_policy = policy;

  if (m_policy == CreationPolicy::GLOBAL) {
    assert(m_callbacks);
    m_cb = std::make_unique<detail::Callback>(m_callbacks,
                                              [](detail::CallbackAction) {});
    m_callbacks->call(make_remote_handle);
    m_cb->operator()(
        CallbackAction{Construct{id(), name, detail::pack(params)}});
  }

  this->do_construct(params);
}

void ObjectHandle::set_parameter(const std::string &name,
                                 const Variant &value) {
  using detail::CallbackAction;
  using SetParameter = CallbackAction::SetParameter;

  if (m_policy == CreationPolicy::GLOBAL) {
    m_cb->operator()(CallbackAction{SetParameter{name, detail::pack(value)}});
  }

  this->do_set_parameter(name, value);
}

Variant ObjectHandle::call_method(const std::string &name,
                                  const VariantMap &params) {
  using detail::CallbackAction;
  using CallMethod = CallbackAction::CallMethod;

  if (m_policy == CreationPolicy::GLOBAL) {
    m_cb->operator()(CallbackAction{CallMethod{name, detail::pack(params)}});
  }

  return this->do_call_method(name, params);
}

void ObjectHandle::initialize(Communication::MpiCallbacks &cb) {
  m_callbacks = &cb;
}
} /* namespace ScriptInterface */
