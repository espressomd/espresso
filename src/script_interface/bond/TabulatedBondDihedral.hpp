#ifndef SCRIPT_INTERFACE_BONDS_TABULATED_BOND_DIHEDRAL_HPP
#define SCRIPT_INTERFACE_BONDS_TABULATED_BOND_DIHEDRAL_HPP

#include "Bond.hpp"

#include "core/bond/TabulatedBondDihedral.hpp"
#include "core/TabulatedPotential.hpp"
#include "get_value.hpp"

#include <utility>
  
namespace ScriptInterface {
namespace Bond {
class TabulatedBondDihedral : public AutoParameters<Bond> {
  using CoreBond = ::Bond::TabulatedBondDihedral;
  std::shared_ptr<CoreBond> m_bond;

public:
  TabulatedBondDihedral() {
    add_parameters(
		   {{"min",
			 [this](Variant const &v) { m_bond->min() = get_value<double>(v); },
			 [this]() { return m_bond->min(); }},
		       {"max",
			   [this](Variant const &v) { m_bond->max() = get_value<double>(v); },
			   [this]() { return m_bond->max(); }},
			 {"energy",
			     [this](Variant const &v) { m_bond->energy() = get_value<std::vector
										     <double> >
			       (v); },
			     [this]() { return m_bond->energy(); }},
			   {"force",
			       [this](Variant const &v) { m_bond->force() = get_value<std::vector
										      <double> >
				 (v); },
			       [this]() { return m_bond->force(); }}
		      
		   });
  }

  std::shared_ptr<::Bond::Bond> bond() const {
    return std::static_pointer_cast<::Bond::Bond>(m_bond);
  }

  void construct(const VariantMap &params) {
    m_bond = make_shared_from_args<CoreBond, double, double, std::vector<double>,
				   std::vector<double> >
      (params, "min", "max", "energy", "force");
  }
  
};
}
}

#endif
