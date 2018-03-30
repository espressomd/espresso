#ifndef SCRIPT_INTERFACE_BONDS_ANGLE_HARMONIC_HPP
#define SCRIPT_INTERFACE_BONDS_ANGLE_HARMONIC_HPP

#include "Bond.hpp"

#include "core/bond/AngleHarmonic.hpp"
#include "get_value.hpp"

#include <utility>
  
namespace ScriptInterface {
namespace Bond {
class AngleHarmonic : public AutoParameters<Bond> {
  using CoreBond = ::Bond::AngleHarmonic;
  std::shared_ptr<CoreBond> m_bond;

public:
  AngleHarmonic() {
    add_parameters(
		   {{"bend",
			 [this](Variant const &v) { m_bond->bend() = get_value<double>(v); },
			 [this]() { return m_bond->bend(); }},
		       {"phi0",
			   [this](Variant const &v) { m_bond->phi0() = get_value<double>(v); },
			   [this]() { return m_bond->phi0(); }}
		   });
  }

  std::shared_ptr<::Bond::Bond> bond() const {
    return std::static_pointer_cast<::Bond::Bond>(m_bond);
  }

  void construct(const VariantMap &params) {
    m_bond = make_shared_from_args<CoreBond, double, double>
      (params, "bend", "phi0");
  }
  
};
}
}

#endif
