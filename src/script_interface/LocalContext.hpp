#ifndef SCRIPT_INTERFACE_LOCAL_CONTEXT_HPP
#define SCRIPT_INTERFACE_LOCAL_CONTEXT_HPP

#include "Context.hpp"
#include "ObjectHandle.hpp"
#include "ObjectState.hpp"

#include <utils/Factory.hpp>

namespace ScriptInterface {
class LocalContext : public Context {
  Utils::Factory<ObjectHandle> m_factory;

public:
  explicit LocalContext(Utils::Factory<ObjectHandle> factory)
      : m_factory(std::move(factory)) {}

  const Utils::Factory<ObjectHandle> &factory() const { return m_factory; }

  void nofity_call_method(const ObjectHandle *, std::string const &,
                          VariantMap const &) override {}
  void notify_set_parameter(const ObjectHandle *, std::string const &,
                            Variant const &) override {}
  void nofity_delete_handle(const ObjectHandle *) override {}

  std::shared_ptr<ObjectHandle>
  make_shared(std::string const &name, const VariantMap &parameters) override {
    auto sp = m_factory.make(name);
    set_manager(sp.get());
    set_name(sp.get(), m_factory.stable_name(name));

    sp->construct(parameters);

    return sp;
  }
};
} // namespace ScriptInterface

#endif
