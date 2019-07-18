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
    m_meta.emplace(so.get(),
                   Meta{CreationPolicy::LOCAL, m_factory.stable_name(name)});
    so->m_manager = shared_from_this();

    so->do_construct(unpack(parameters, m_local_objects));

    m_local_objects.emplace(std::make_pair(id, std::move(so)));
  } catch (std::runtime_error const &) {
  }
}

void ObjectManager::remote_make_handle(ObjectId id, const std::string &name,
                                       const VariantMap &parameters) {
  cb_make_handle(id, name, pack(parameters));
}

void ObjectManager::remote_delete_handle(const ObjectHandle *o) {
  if (policy(o) == CreationPolicy::GLOBAL) {
    cb_delete_handle(object_id(o));
  }
}

void ObjectManager::set_parameter(ObjectId id, std::string const &name,
                                  PackedVariant const &value) {
  try {
    m_local_objects.at(id)->do_set_parameter(name,
                                             unpack(value, m_local_objects));
  } catch (std::runtime_error const &) {
  }
}

void ObjectManager::remote_set_parameter(const ObjectHandle *o,
                                         std::string const &name,
                                         Variant const &value) {
  if (policy(o) == CreationPolicy::GLOBAL) {
    cb_set_parameter(object_id(o), name, pack(value));
  }
}

void ObjectManager::call_method(ObjectId id, std::string const &name,
                                PackedMap const &arguments) {
  try {
    m_local_objects.at(id)->do_call_method(name,
                                           unpack(arguments, m_local_objects));
  } catch (std::runtime_error const &) {
  }
}

void ObjectManager::remote_call_method(const ObjectHandle *o,
                                       std::string const &name,
                                       VariantMap const &arguments) {
  if (policy(o) == CreationPolicy::GLOBAL) {
    cb_call_method(object_id(o), name, pack(arguments));
  }
}

std::shared_ptr<ObjectHandle>
ObjectManager::make_shared(std::string const &name, CreationPolicy policy,
                           const VariantMap &parameters) {
  auto sp = m_factory.make(name);

  auto const id = object_id(sp.get());

  sp->m_manager = shared_from_this();

  m_meta.emplace(sp.get(), Meta{policy, m_factory.stable_name(name)});

  if (policy == CreationPolicy::GLOBAL) {
    remote_make_handle(id, name, parameters);
  }

  sp->do_construct(parameters);

  return sp;
}

std::string ObjectManager::serialize(const ObjectRef &o) const {
  ObjectState state{std::string{name(o.get())},
                    policy(o.get()),
                    {},
                    {},
                    o->get_internal_state()};

  auto const params = o->get_parameters();
  state.params.resize(params.size());

  PackVisitor v;

  /* Pack parameters and keep track of ObjectRef parameters */
  boost::transform(params, state.params.begin(),
                   [&v](auto const &kv) -> PackedMap::value_type {
                     return {kv.first, boost::apply_visitor(v, kv.second)};
                   });

  /* Packed Object parameters */
  state.objects.resize(v.objects().size());
  boost::transform(v.objects(), state.objects.begin(), [this](auto const &kv) {
    return std::make_pair(kv.first, serialize(kv.second));
  });

  return Utils::pack(state);
}

std::shared_ptr<ObjectHandle>
ObjectManager::unserialize(std::string const &state_) {
  auto state = Utils::unpack<ObjectState>(state_);

  std::unordered_map<ObjectId, ObjectRef> objects;
  boost::transform(state.objects, std::inserter(objects, objects.end()),
                   [this](auto const &kv) {
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
} // namespace ScriptInterface