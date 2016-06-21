#ifndef SCRIPT_INTERFACE_CONSTRAINTS_CHARGED_ROD_HPP
#define SCRIPT_INTERFACE_CONSTRAINTS_CHARGED_ROD_HPP

#include "ScriptInterface.hpp"
#include "core/constraints/ChargedRod.hpp"

namespace ScriptInterface {
namespace Constraints {

class ChargedRod : public ScriptInterfaceBase {
public:
  /** Parsing stuff */
  VariantMap get_parameters() const override;
  ParameterMap all_parameters() const override;
  void set_parameter(const std::string &name, const Variant &value) override;
  const std::string name() const override { return "Constraints::ChargedRod"; }

private:
  ::Constraints::ChargedRod m_rod;
};
}
}

#endif
