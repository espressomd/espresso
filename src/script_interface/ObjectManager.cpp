#include "ObjectManager.hpp"
#include "ObjectHandle.hpp"
#include "ObjectState.hpp"
#include "PackedVariant.hpp"

#include <utils/serialization/pack.hpp>

namespace ScriptInterface {
void ObjectManager::make_handle(ObjectId id, const std::string &name,
                                const PackedMap &parameters) {
  try {
    m_local_objects[id] =
        make_shared(name, CreationPolicy::LOCAL,
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
  auto sp = m_factory.make(name);

  sp->m_manager = this;
  sp->m_name = name;
  sp->m_policy = policy;

  if (policy == CreationPolicy::GLOBAL) {
    remote_make_handle(object_id(sp.get()), name, parameters);
  }

  sp->do_construct(parameters);

  return sp;
}

std::string ObjectManager::serialize(const ObjectRef &o) const {
  ObjectState state{o->name(), o->policy(), {}, {}, o->get_internal_state()};

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