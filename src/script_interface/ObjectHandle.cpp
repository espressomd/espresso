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
#include <boost/serialization/variant.hpp>

#include <sstream>

namespace ScriptInterface {
namespace detail {
struct CallbackAction {
  struct Construct {
    std::string name;
    VariantMap parameters;

    template <class Archive> void serialize(Archive &ar, long int) {
      ar &name &parameters;
    }
  };
  struct SetParameter {
    std::string name;
    Variant value;

    template <class Archive> void serialize(Archive &ar, long int) {
      ar &name &value;
    }
  };
  struct CallMethod {
    std::string name;
    VariantMap arguments;

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
      o = ObjectHandle::make_shared(
          ctor.name, ObjectHandle::CreationPolicy::LOCAL, ctor.parameters);
    }
    void operator()(const CallbackAction::SetParameter &param) const {
      assert(o), o->set_parameter(param.name, param.value);
    }
    void operator()(const CallbackAction::CallMethod &method) const {
      assert(o);
      (void)o->call_method(method.name, method.arguments);
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

  /* Id of the newly created instance */
  const auto id = sp->id();

  /* Now get a reference to the corresponding weak_ptr in ObjectId and update
     it with our shared ptr, so that everybody uses the same ref count.
  */
  sp->get_instance(id) = sp;

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

  switch (policy) {
  case CreationPolicy::LOCAL:
    break;
  case CreationPolicy::GLOBAL:
    assert(m_callbacks);
    m_cb_ = std::make_unique<detail::Callback>(m_callbacks,
                                               [](detail::CallbackAction) {});
    m_callbacks->call(make_remote_handle);
    m_cb_->operator()(CallbackAction{Construct{name, params}});
  }

  this->do_construct(params);
}

void ObjectHandle::set_parameter(const std::string &name,
                                 const Variant &value) {
  this->do_set_parameter(name, value);
}

Variant ObjectHandle::call_method(const std::string &name,
                                  const VariantMap &params) {
  return this->do_call_method(name, params);
}

void ObjectHandle::initialize(Communication::MpiCallbacks &cb) {
  m_callbacks = &cb;
}
} /* namespace ScriptInterface */
