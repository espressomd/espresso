#include "ObjectManager.hpp"
#include "ObjectHandle.hpp"
#include "PackedVariant.hpp"

#include <utils/serialization/pack.hpp>

namespace ScriptInterface {
void ObjectManager::make_handle(ObjectId id, const std::string &name,
                                const PackedMap &parameters) {
  try {
    m_local_objects[id] =
        ObjectHandle::make_shared(name, CreationPolicy::LOCAL,
                                  unpack(parameters, m_local_objects));
  } catch (std::runtime_error const &) {
  }
}

void ObjectManager::remote_make_handle(ObjectId id, const std::string &name,
                                       const VariantMap &parameters) {
  cb_make_handle(id, name, pack(parameters));
}

void ObjectManager::set_parameter(ObjectId id, std::string const &name,
                                  PackedVariant const &value) {
  try {
    m_local_objects.at(id)->set_parameter(name, unpack(value, m_local_objects));
  } catch (std::runtime_error const &) {
  }
}

void ObjectManager::remote_set_parameter(ObjectId id, std::string const &name,
                                         Variant const &value) {
  cb_set_parameter(id, name, pack(value));
}

void ObjectManager::call_method(ObjectId id, std::string const &name,
                                PackedMap const &arguments) {
  try {
    m_local_objects.at(id)->call_method(name, unpack(arguments, m_local_objects));
  } catch (std::runtime_error const &) {
  }
}

void ObjectManager::remote_call_method(ObjectId id, std::string const &name,
                                       VariantMap const &arguments) {
  cb_call_method(id, name, pack(arguments));
}

std::shared_ptr<ObjectHandle>
ObjectManager::make_shared(std::string const &name,
                           CreationPolicy policy,
                           const VariantMap &parameters) {
  auto sp = factory.make(name);

  sp->m_manager = this;
  sp->m_name = name;
  sp->m_policy = policy;

  if (sp->m_policy == CreationPolicy::GLOBAL) {
    remote_make_handle(object_id(sp.get()), name, parameters);
  }

  sp->do_construct(parameters);

  return sp;
}
} // namespace ScriptInterface