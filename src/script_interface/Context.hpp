#ifndef ESPRESSO_CONTEXT_HPP
#define ESPRESSO_CONTEXT_HPP

#include "ObjectHandle.hpp"
#include "ObjectState.hpp"
#include "Variant.hpp"

#include <boost/utility/string_ref.hpp>
#include <utils/serialization/pack.hpp>

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
  virtual void nofity_call_method(const ObjectHandle *self,
                                  std::string const &name,
                                  VariantMap const &arguments) = 0;

  /**
   * @brief Set a parameter on remote instances
   *
   * @param o Internal identifier of the instance to be modified
   * @param name Name of the parameter to change
   * @param value Value to set it to
   */
  virtual void notify_set_parameter(const ObjectHandle *self,
                                    std::string const &name,
                                    Variant const &value) = 0;

  /**
   * @brief Delete remote instances
   *
   * @param o Internal identified of the instance
   */
  virtual void nofity_delete_handle(const ObjectHandle *self) = 0;

  /**
   * @brief Get a new reference counted instance of a script interface by
   * name.
   *
   */
  virtual std::shared_ptr<ObjectHandle>
  make_shared(std::string const &name, const VariantMap &parameters) = 0;

  /**
   * @brief String representation of the state of an object.
   */
  std::string serialize(const ObjectHandle *o) const {
    ObjectState state;

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
    boost::transform(
        v.objects(), state.objects.begin(), [this](auto const &kv) {
          return std::make_pair(kv.first, serialize(kv.second.get()));
        });

    state.name = name(o).to_string();
    state.internal_state = o->get_internal_state();

    return Utils::pack(state);
  }

  /**
   * @brief Make object from serialized state.
   */
  ObjectRef deserialize(std::string const &state_) {
    auto const state = Utils::unpack<ObjectState>(state_);

    std::unordered_map<ObjectId, ObjectRef> objects;
    boost::transform(state.objects, std::inserter(objects, objects.end()),
                     [this](auto const &kv) {
                       return std::make_pair(kv.first,
                                             deserialize(kv.second));
                     });

    VariantMap params;
    for (auto const &kv : state.params) {
      params[kv.first] =
          boost::apply_visitor(UnpackVisitor(objects), kv.second);
    }

    auto o = make_shared(state.name, {});
    o->construct(params);
    o->set_internal_state(state.internal_state);

    return o;
  }

protected:
  void set_manager(ObjectHandle *o) { o->m_manager = this->shared_from_this(); }

  void set_name(ObjectHandle *o, boost::string_ref name) const {
    o->m_name = name;
  }

public:
  boost::string_ref name(const ObjectHandle *o) const { return o->name(); }

public:
  virtual ~Context() = default;
};
} // namespace ScriptInterface
#endif // ESPRESSO_CONTEXT_HPP
