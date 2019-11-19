#include "ObjectManager.hpp"

#include "GlobalContext.hpp"
#include "LocalContext.hpp"

namespace ScriptInterface {
std::shared_ptr<ObjectHandle>
ObjectManager::make_shared(CreationPolicy policy, std::string const &name,
                           const VariantMap &parameters) {
  return context(policy)->make_shared(name, parameters);
}

std::shared_ptr<ObjectHandle>
ObjectManager::deserialize(std::string const &state_) {
  auto const state =
      Utils::unpack<std::pair<CreationPolicy, std::string>>(state_);

  return context(state.first)->deserialize(state.second);
}

std::string ObjectManager::serialize(const ObjectHandle *o) const {
  auto ctx = o->manager();
  return assert(ctx),
         Utils::pack(std::make_pair(policy(ctx), ctx->serialize(o)));
}

ObjectManager::ObjectManager(Communication::MpiCallbacks &callbacks,
                             const Utils::Factory<ObjectHandle> &factory)
    : m_local_context(std::make_shared<LocalContext>(factory)),
      m_global_context(std::make_shared<GlobalContext>(callbacks, factory)) {}
} // namespace ScriptInterface
