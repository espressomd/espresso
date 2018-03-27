#ifndef SCRIPT_INTERFACE_BONDS_BONDED_COULOMB_HPP
#define SCRIPT_INTERFACE_BONDS_BONDED_COULOMB_HPP

#include "Bond.hpp"

#include "core/bond/BondedCoulomb.hpp"
#include "get_value.hpp"

#include <utility>
  
namespace ScriptInterface {
namespace Bond {
class BondedCoulomb : public AutoParameters<Bond> {
  using CoreBond = ::Bond::BondedCoulomb;
  std::shared_ptr<CoreBond> m_bond;

public:
 BondedCoulomb () {

   add_parameters(
		  {{"prefactor",
			[this](Variant const &v) { m_bond->prefactor() = get_value<double>(v); },
			[this]() { return m_bond->prefactor(); }},
		      });
 }
      

  std::shared_ptr<::Bond::Bond> bond() const {
    return std::static_pointer_cast<::Bond::Bond>(m_bond);
  }

  void construct(const VariantMap &params) {
    m_bond = make_shared_from_args<CoreBond, double>(params, "prefactor");
  }
 };
}
}

#endif
