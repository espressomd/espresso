#ifndef SCRIPT_INTERFACE_BONDS_BONDED_COULOMB_P3MSR_HPP
#define SCRIPT_INTERFACE_BONDS_BONDED_COULOMB_P3MSR_HPP

#include "Bond.hpp"

#include "core/bond/BondedCoulombP3MSR.hpp"
#include "get_value.hpp"

#include <utility>
  
namespace ScriptInterface {
namespace Bond {
class BondedCoulombP3MSR : public AutoParameters<Bond> {
  using CoreBond = ::Bond::BondedCoulombP3MSR;
  std::shared_ptr<CoreBond> m_bond;

public:
 BondedCoulombP3MSR () {

   add_parameters(
		  {{"q1q2",
			[this](Variant const &v) { m_bond->q1q2() = get_value<double>(v); },
			[this]() { return m_bond->q1q2(); }},
		      });
 }
      

  std::shared_ptr<::Bond::Bond> bond() const {
    return std::static_pointer_cast<::Bond::Bond>(m_bond);
  }

  void construct(const VariantMap &params) {
    m_bond = make_shared_from_args<CoreBond, double>(params, "q1q2");
  }
 };
}
}

#endif
