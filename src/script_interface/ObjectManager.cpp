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
                             const Utils::Factory<ObjectHandle> &factory) {
  auto local_context = std::make_shared<LocalContext>(factory);

  /* If there is only one node, we can treat all objects as local, and thus
   * never invoke any callback. */
  m_global_context =
      (callbacks.comm().size() > 1)
          ? std::make_shared<GlobalContext>(callbacks, local_context)
          : std::static_pointer_cast<Context>(local_context);

  m_local_context  = std::move(local_context);
}
} // namespace ScriptInterface
