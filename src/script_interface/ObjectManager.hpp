#ifndef ESPRESSO_SCRIPT_INTERFACE_OBJECTMANAGER_HPP
#define ESPRESSO_SCRIPT_INTERFACE_OBJECTMANAGER_HPP

#include "MpiCallbacks.hpp"
#include "ObjectHandle.hpp"
#include "PackedVariant.hpp"

#include <utils/Factory.hpp>

#include <boost/serialization/utility.hpp>

namespace ScriptInterface {
extern Utils::Factory<ObjectHandle> factory;
} /* namespace ScriptInterface */

namespace ScriptInterface {
class ObjectManager {
  using ObjectId = std::size_t;

  /* Instances on this node that are managed by the
   * head node. */
  std::unordered_map<ObjectId, ObjectRef> m_local_objects;

public:
  auto const &local_objects() const { return m_local_objects; }

private:
  Communication::CallbackHandle<ObjectId, const std::string &,
                                const PackedMap &>
      cb_make_handle;
  Communication::CallbackHandle<ObjectId, const std::string &,
                                const PackedVariant &>
      cb_set_parameter;
  Communication::CallbackHandle<ObjectId, std::string const &,
                                PackedMap const &>
      cb_call_method;
  Communication::CallbackHandle<ObjectId> cb_delete_handle;

public:
  explicit ObjectManager(Communication::MpiCallbacks *callbacks)
      : cb_make_handle(callbacks,
                       [this](ObjectId id, const std::string &name,
                              const PackedMap &parameters) {
                         make_handle(id, name, parameters);
                       }),
        cb_set_parameter(callbacks,
                         [this](ObjectId id, std::string const &name,
                                PackedVariant const &value) {
                           set_parameter(id, name, value);
                         }),
        cb_call_method(callbacks,
                       [this](ObjectId id, std::string const &name,
                              PackedMap const &arguments) {
                         call_method(id, name, arguments);
                       }),
        cb_delete_handle(callbacks,
                         [this](ObjectId id) { delete_handle(id); }) {}

private:
  /**
   * @brief Callback for @function remote_make_handle
   */
  void make_handle(ObjectId id, const std::string &name,
                   const PackedMap &parameters);

public:
  /**
   * @brief Create remote instances
   *
   * @param id Internal identifier of the instance
   * @param name Class name
   * @param parameters Constructor parameters.
   */
  void remote_make_handle(ObjectId id, const std::string &name,
                          const VariantMap &parameters);

private:
  /**
   * @brief Callback for @function remote_set_parameter
   */
  void set_parameter(ObjectId id, std::string const &name,
                     PackedVariant const &value);

public:
  /**
   * @brief Set a parameter on remote instances
   *
   * @param id Internal identifier of the instance to be modified
   * @param name Name of the parameter to change
   * @param value Value to set it to
   */
  void remote_set_parameter(ObjectId id, std::string const &name,
                            Variant const &value);

private:
  /**
   * @brief Callback for @function remote_call_method
   */
  void call_method(ObjectId id, std::string const &name,
                   PackedMap const &arguments);

public:
  /**
   * @brief Call method on remote instances
   *
   * @param id Internal identified of the instance
   * @param name Name of the method to call
   * @param arguments Arguments to the call
   */
  void remote_call_method(ObjectId id, std::string const &name,
                          VariantMap const &arguments);

private:
  /**
   * @brief Callback for @function remote_delete_handle
   */
  void delete_handle(ObjectId id) { m_local_objects.erase(id); }

public:
  /**
   * @brief Delete remote instances
   *
   * @param id Internal identified of the instance
   */
  void remote_delete_handle(ObjectId id) { cb_delete_handle(id); }

  template <typename T> static void register_new(std::string const &name) {
    static_assert(std::is_base_of<ObjectHandle, T>::value, "");

    /* Register with the factory */
    factory.register_new<T>(name);
  }

  /**
   * @brief Get a new reference counted instance of a script interface by
   * name.
   *
   */
  std::shared_ptr<ObjectHandle> make_shared(std::string const &name,
                                            CreationPolicy policy,
                                            const VariantMap &parameters = {});

  /**
   * @brief Returns a binary representation of the state of a object.
   */
  std::string serialize(const ObjectRef &o) const;

  /**
   * @brief Creates a new instance from a binary state,
   *        as returned by @function serialize.
   */
  std::shared_ptr<ObjectHandle> unserialize(std::string const &state_);
};
} // namespace ScriptInterface

#endif
