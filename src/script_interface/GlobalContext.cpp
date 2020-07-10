#include "GlobalContext.hpp"
#include "ObjectHandle.hpp"
#include "PackedVariant.hpp"

#include <utils/serialization/pack.hpp>

namespace ScriptInterface {
void GlobalContext::make_handle(ObjectId id, const std::string &name,
                                const PackedMap &parameters) {
  try {
    ObjectRef so = m_node_local_context->make_shared(
        name, unpack(parameters, m_local_objects));

    m_local_objects.emplace(std::make_pair(id, std::move(so)));
  } catch (std::runtime_error const &) {
  }
}

void GlobalContext::remote_make_handle(ObjectId id, const std::string &name,
                                       const VariantMap &parameters) {
  cb_make_handle(id, name, pack(parameters));
}

void GlobalContext::set_parameter(ObjectId id, std::string const &name,
                                  PackedVariant const &value) {
  try {
    m_local_objects.at(id)->set_parameter(name, unpack(value, m_local_objects));
  } catch (std::runtime_error const &) {
  }
}

void GlobalContext::notify_set_parameter(const ObjectHandle *o,
                                         std::string const &name,
                                         Variant const &value) {
  cb_set_parameter(object_id(o), name, pack(value));
}

void GlobalContext::call_method(ObjectId id, std::string const &name,
                                PackedMap const &arguments) {
  try {
    m_local_objects.at(id)->call_method(name,
                                        unpack(arguments, m_local_objects));
  } catch (std::runtime_error const &) {
  }
}

void GlobalContext::notify_call_method(const ObjectHandle *o,
                                       std::string const &name,
                                       VariantMap const &arguments) {
  cb_call_method(object_id(o), name, pack(arguments));
}

std::shared_ptr<ObjectHandle>
GlobalContext::make_shared(std::string const &name,
                           const VariantMap &parameters) {
  std::unique_ptr<ObjectHandle> sp = m_node_local_context->factory().make(name);
  set_manager(sp.get());
  set_name(sp.get(), m_node_local_context->factory().stable_name(name));

  auto const id = object_id(sp.get());
  remote_make_handle(id, name, parameters);

  sp->construct(parameters);

  return {sp.release(),
          /* Custom deleter, we keep the corresponding global context,
           * as well as the original deleter for the object. */
          [global_context = this, deleter = sp.get_deleter()](ObjectHandle *o) {
            /* Tell the other nodes before invoking the destructor, this is
             * required
             * to have synchronous destructors, which is needed by some client
             * code. */
            global_context->cb_delete_handle(object_id(o));

            /* Locally destroy the object. */
            deleter(o);
          }};
}
} // namespace ScriptInterface