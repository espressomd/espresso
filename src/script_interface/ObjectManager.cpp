#include "ObjectManager.hpp"
#include "ObjectHandle.hpp"
#include "PackedVariant.hpp"

#include <utils/serialization/pack.hpp>

#include <boost/serialization/utility.hpp>

namespace ScriptInterface {
void ObjectManager::make_handle(ObjectId id, const std::string &name,
                 const PackedMap &parameters) {
  try {
    local_objects[id] =
        ObjectHandle::make_shared(name, ObjectHandle::CreationPolicy::LOCAL,
                                  unpack(parameters, local_objects));
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
    local_objects.at(id)->set_parameter(name, unpack(value, local_objects));
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
    local_objects.at(id)->call_method(name, unpack(arguments, local_objects));
  } catch (std::runtime_error const &) {
  }
}

void ObjectManager::remote_call_method(ObjectId id, std::string const &name,
                        VariantMap const &arguments) {
  cb_call_method(id, name, pack(arguments));
}

}