#ifndef ESPRESSO_OBJECTMANAGER_HPP
#define ESPRESSO_OBJECTMANAGER_HPP

#include "Context.hpp"
#include "Variant.hpp"

#include "MpiCallbacks.hpp"

#include <utils/Factory.hpp>

#include <memory>
#include <stdexcept>

namespace ScriptInterface {
class ObjectManager {
  std::shared_ptr<Context> m_local_context;
  std::shared_ptr<Context> m_global_context;

public:
  enum class CreationPolicy { LOCAL, GLOBAL };

  ObjectManager(Communication::MpiCallbacks &callbacks,
                const Utils::Factory<ObjectHandle> &factory);

  /**
   * @brief Get a new reference counted instance of a script interface by
   * name.
   */
  std::shared_ptr<ObjectHandle> make_shared(CreationPolicy policy,
                                            std::string const &name,
                                            const VariantMap &parameters);

  /**
   * @brief Get a new reference counted instance of a script interface from
   *         a serialized state.
   */
  std::shared_ptr<ObjectHandle> deserialize(std::string const &state_);

  /**
   * @brief Serialize a script interface object into a binary representation.
   */
  std::string serialize(const ObjectHandle *o) const;

private:
  /**
   * @brief Map policy to context.
   *
   * Inverse of policy.
   */
  Context *context(CreationPolicy policy) const {
    switch (policy) {
    case CreationPolicy::LOCAL:
      return assert(m_local_context), m_local_context.get();
    case CreationPolicy::GLOBAL:
      return assert(m_global_context), m_global_context.get();
    default:
      throw std::runtime_error("Unknown context type.");
    }
  }

  /**
   * @brief Map context to policy.
   *
   * Inverse of context.
   */
  CreationPolicy policy(Context *c) const {
    if (c == m_local_context.get()) {
      return CreationPolicy::LOCAL;
    } else if (c == m_global_context.get()) {
      return CreationPolicy::GLOBAL;
    }

    throw std::runtime_error("Invalid context.");
  }
};
} // namespace ScriptInterface

#endif // ESPRESSO_OBJECTMANAGER_HPP
