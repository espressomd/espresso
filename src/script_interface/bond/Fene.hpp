#ifndef SCRIPT_INTERFACE_BONDS_FENE_HPP
#define SCRIPT_INTERFACE_BONDS_FENE_HPP

#include "Bond.hpp"

#include "core/bond/Fene.hpp"
#include "get_value.hpp"

#include <utility>
  
namespace ScriptInterface {
namespace Bond {
class Fene : public AutoParameters<Bond> {
  using CoreBond = ::Bond::Fene;
  std::shared_ptr<CoreBond> m_bond;

public:
  Fene() {
    add_parameters(
        {{"r0",
          [this](Variant const &v) { m_bond->r0() = get_value<double>(v); },
          [this]() { return m_bond->r0(); }},
{"k",
          [this](Variant const &v) { m_bond->k() = get_value<double>(v); },
          [this]() { return m_bond->k(); }},
{"dr_max",
    [this](Variant const &v) { m_bond->set_dr_max(  get_value<double>(v)) ; },
          [this]() { return m_bond->dr_max(); }}
});
  }

  std::shared_ptr<::Bond::Bond> bond() const {
    return std::static_pointer_cast<::Bond::Bond>(m_bond);
  }

  void construct(const VariantMap &params) {
    m_bond = make_shared_from_args<CoreBond, double, double, double>(
        params, "r0", "dr_max", "k");
  }
};
}
}

#endif
