#ifndef ESPRESSO_CONTEXT_HPP
#define ESPRESSO_CONTEXT_HPP

#include "Variant.hpp"

#include <utils/Factory.hpp>

#include <boost/utility/string_ref.hpp>

#include <memory>

namespace ScriptInterface {
class Context : public std::enable_shared_from_this<Context> {
public:
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
  virtual std::string serialize(const ObjectRef &o) const;

  /**
    * @brief Initialize a bare object from binary state,
    * as returned by @function serialize.
    */
  void
  unserialize(ObjectHandle *o, std::string const &state_) const;

public:
  /**
   * @brief Register new class type for this context.
   *
   * @tparam T class type derived from @class ObjectHandle.
   * @param name Name by which the class should be registered.
   */
  template <typename T>
  std::enable_if_t<std::is_base_of<ObjectHandle, T>::value>
      register_new(std::string const &name) {
    /* Register with the factory */
    m_factory.register_new<T>(name);
  }

  virtual std::shared_ptr<ObjectHandle> make_shared(const ObjectHandle *,
                                            std::string const &name,
                                            const VariantMap &parameters) {};

protected:
  /**
   * @brief Make a bare new object
   *
   * @param name Class name
   * @return New instance
   */
  ObjectRef make_bare(const std::string &name);
public:
  virtual ~Context() = default;

  boost::string_ref name(const ObjectHandle *) const;
private:
  Utils::Factory<ObjectHandle> m_factory;
};
}
#endif //ESPRESSO_CONTEXT_HPP
