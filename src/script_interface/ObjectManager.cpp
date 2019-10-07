#include "ObjectManager.hpp"
#include "ObjectHandle.hpp"
#include "ObjectState.hpp"
#include "PackedVariant.hpp"

#include <utils/serialization/pack.hpp>

namespace ScriptInterface {
void ObjectManager::make_handle(ObjectId id, const std::string &name,
                                const PackedMap &parameters) {
  try {
    ObjectRef so = m_factory.make(name);
    m_policy.emplace(so.get(), CreationPolicy::LOCAL);

    so->construct(unpack(parameters, m_local_objects));

    m_local_objects.emplace(std::make_pair(id, std::move(so)));
  } catch (std::runtime_error const &) {
  }
}

void ObjectManager::remote_make_handle(ObjectId id, const std::string &name,
                                       const VariantMap &parameters) {
  cb_make_handle(id, name, pack(parameters));
}

void ObjectManager::nofity_delete_handle(const ObjectHandle *o) {
  if (policy(o) == CreationPolicy::GLOBAL) {
    cb_delete_handle(object_id(o));
  }
}

void ObjectManager::set_parameter(ObjectId id, std::string const &name,
                                  PackedVariant const &value) {
  try {
    m_local_objects.at(id)->set_parameter(name, unpack(value, m_local_objects));
  } catch (std::runtime_error const &) {
  }
}

void ObjectManager::notify_set_parameter(const ObjectHandle *o,
                                         std::string const &name,
                                         Variant const &value) {
  if (policy(o) == CreationPolicy::GLOBAL) {
    cb_set_parameter(object_id(o), name, pack(value));
  }
}

void ObjectManager::call_method(ObjectId id, std::string const &name,
                                PackedMap const &arguments) {
  try {
    m_local_objects.at(id)->call_method(name,
                                        unpack(arguments, m_local_objects));
  } catch (std::runtime_error const &) {
  }
}

void ObjectManager::nofity_call_method(const ObjectHandle *o,
                                       std::string const &name,
                                       VariantMap const &arguments) {
  if (policy(o) == CreationPolicy::GLOBAL) {
    cb_call_method(object_id(o), name, pack(arguments));
  }
}

std::shared_ptr<ObjectHandle>
ObjectManager::make_shared(const ObjectHandle *, std::string const &name,
                           CreationPolicy policy,
                           const VariantMap &parameters) {
  auto sp = m_factory.make(name);

  auto const id = object_id(sp.get());

  m_policy.emplace(sp.get(), policy);

  if (policy == CreationPolicy::GLOBAL) {
    remote_make_handle(id, name, parameters);
  }

  sp->construct(parameters);

  return sp;
}

std::shared_ptr<ObjectHandle>
ObjectManager::make_shared(const ObjectHandle *self, std::string const &name,
                           const VariantMap &parameters) {
  return make_shared(self, name, policy(self), parameters);
}
} // namespace ScriptInterface