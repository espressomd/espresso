#include "Context.hpp"
#include "ObjectHandle.hpp"
#include "ObjectState.hpp"

#include <utils/serialization/pack.hpp>

#include <boost/serialization/utility.hpp>

namespace ScriptInterface {
ObjectRef Context::make_bare(const std::string &name) {
  auto sp = m_factory.make(name);

  sp->m_manager = shared_from_this();
  sp->m_name = m_factory.stable_name(name);

  return sp;
}

boost::string_ref Context::name(const ObjectHandle *o) const {
  return o->m_name;
}

std::string Context::serialize(const ObjectRef &o) const {
  ObjectState state{
      std::string{name(o.get())}, {}, {}, o->get_internal_state()};

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

void Context::unserialize(ObjectHandle *o, std::string const &state_) const {
  /*
  auto const state = Utils::unpack<ObjectState>(state_);

  std::unordered_map<ObjectId, ObjectRef> objects;
  boost::transform(state.objects, std::inserter(objects, objects.end()),
                   [this, o](auto const &kv) {
                     return std::make_pair(kv.first, unserialize(o, kv.second));
                   });

  VariantMap params;
  for (auto const &kv : state.params) {
    params[kv.first] = boost::apply_visitor(UnpackVisitor(objects), kv.second);
  }

  o->construct(params);
  o->set_internal_state(state.internal_state);
  */
}
}