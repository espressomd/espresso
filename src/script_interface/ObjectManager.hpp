#ifndef ESPRESSO_SCRIPT_INTERFACE_OBJECTMANAGER_HPP
#define ESPRESSO_SCRIPT_INTERFACE_OBJECTMANAGER_HPP

#include "MpiCallbacks.hpp"
#include "ObjectHandle.hpp"
#include "PackedVariant.hpp"

#include <utils/Factory.hpp>

#include <boost/serialization/utility.hpp>

namespace ScriptInterface {

class Context {
  /**
   * @brief Call method on remote instances
   *
   * @param o Internal identified of the instance
   * @param name Name of the method to call
   * @param arguments Arguments to the call
   */
  virtual void nofity_call_method(const ObjectHandle *o,
                                  std::string const &name,
                                  VariantMap const &arguments) = 0;

  /**
   * @brief Set a parameter on remote instances
   *
   * @param o Internal identifier of the instance to be modified
   * @param name Name of the parameter to change
   * @param value Value to set it to
   */
  virtual void notify_set_parameter(const ObjectHandle *o,
                                    std::string const &name,
                                    Variant const &value) = 0;

/**
   * @brief Delete remote instances
   *
   * @param o Internal identified of the instance
   */
  virtual void nofity_delete_handle(const ObjectHandle *o) = 0;

  /**
   * @brief Returns a binary representation of the state of a object.
   */
  virtual std::string serialize(const ObjectRef &o) const = 0;

  /**
   * @brief Creates a new instance from a binary state,
   *        as returned by @function serialize.
   */
  virtual std::shared_ptr<ObjectHandle> unserialize(std::string const &state_) = 0;
  
public:
  template <typename T> void register_new(std::string const &name) {
    static_assert(std::is_base_of<ObjectHandle, T>::value, "");

    /* Register with the factory */
    m_factory.register_new<T>(name);
  }

  virtual ~Context() = default;
  auto const& factory() const {
    return m_factory;
  }

private:
  Utils::Factory<ObjectHandle> m_factory;
};

class ObjectManager : public Context,
                      public std::enable_shared_from_this<ObjectManager> {
  using ObjectId = std::size_t;

  struct Meta {
    CreationPolicy policy;
    boost::string_ref class_name;
  };

  /* Instances on this node that are managed by the
   * head node. */
  std::unordered_map<ObjectId, ObjectRef> m_local_objects;
  std::unordered_map<const ObjectHandle *, Meta> m_meta;

  const Meta &meta(const ObjectHandle *o) const { return m_meta.at(o); }

  CreationPolicy policy(const ObjectHandle *o) const { return meta(o).policy; }

public:
  boost::string_ref name(const ObjectHandle *o) const {
    return meta(o).class_name;
  }

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
  explicit ObjectManager(Communication::MpiCallbacks &callbacks)
      : cb_make_handle(&callbacks,
                       [this](ObjectId id, const std::string &name,
                              const PackedMap &parameters) {
                         make_handle(id, name, parameters);
                       }),
        cb_set_parameter(&callbacks,
                         [this](ObjectId id, std::string const &name,
                                PackedVariant const &value) {
                           set_parameter(id, name, value);
                         }),
        cb_call_method(&callbacks,
                       [this](ObjectId id, std::string const &name,
                              PackedMap const &arguments) {
                         call_method(id, name, arguments);
                       }),
        cb_delete_handle(&callbacks,
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
  void notify_set_parameter(const ObjectHandle *o, std::string const &name,
                            Variant const &value) override;

private:
  /**
   * @brief Callback for @function remote_call_method
   */
  void call_method(ObjectId id, std::string const &name,
                   PackedMap const &arguments);

public:
  void nofity_call_method(const ObjectHandle *o, std::string const &name,
                          VariantMap const &arguments) override;

private:
  /**
   * @brief Callback for @function remote_delete_handle
   */
  void delete_handle(ObjectId id) { m_local_objects.erase(id); }

public:
  /**
   * @brief Delete remote instances
   *
   * @param o Internal identified of the instance
   */
  void nofity_delete_handle(const ObjectHandle *o) override;

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
  std::string serialize(const ObjectRef &o) const override;

  /**
   * @brief Creates a new instance from a binary state,
   *        as returned by @function serialize.
   */
  std::shared_ptr<ObjectHandle> unserialize(std::string const &state_) override;
};
} // namespace ScriptInterface

#endif
