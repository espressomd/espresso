#ifndef SCRIPT_INTERFACE_BONDS_OVERLAP_BOND_DIHEDRAL_HPP
#define SCRIPT_INTERFACE_BONDS_OVERLAP_BOND_DIHEDRAL_HPP

#include "Bond.hpp"

#include "core/bond/OverlapBondDihedral.hpp"
#include "core/bond/OverlappedBondedInteraction.hpp"
#include "get_value.hpp"

#include <utility>
  
namespace ScriptInterface {
namespace Bond {
class OverlapBondDihedral : public AutoParameters<Bond> {
  using CoreBond = ::Bond::OverlapBondDihedral;
  std::shared_ptr<CoreBond> m_bond;

public:
  OverlapBondDihedral() {
    add_parameters(
		   {{"filename",
			 [this](Variant const &v) {m_bond->filename() = get_value<std::string>(v);},
			 [this]() { return m_bond->filename(); }},
		       {"maxval",
			   [this](Variant const &v) { m_bond->maxval() = get_value<double>(v); },
			   [this]() { return m_bond->maxval(); }},		       
			 {"noverlaps",
			     [this](Variant const &v) { m_bond->noverlaps() = get_value<int>(v); },
			     [this]() { return m_bond->noverlaps(); }},
			   {"para_a",
			       [this](Variant const &v) {m_bond->para_a() = get_value<std::vector
										      <double>>(v);},
			       [this]() { return m_bond->para_a(); }},
			     {"para_b",
				 [this](Variant const &v) {m_bond->para_b() = get_value<std::vector
											  <double>>
				   (v);},
				 [this]() { return m_bond->para_b(); }},
			       {"para_c",
				   [this](Variant const &v) {m_bond->para_c() = get_value<std::vector
											  <double>>
				     (v); },
				   [this]() { return m_bond->para_c(); }}
		      
		   });
  }

  std::shared_ptr<::Bond::Bond> bond() const {
    return std::static_pointer_cast<::Bond::Bond>(m_bond);
  }

  void construct(const VariantMap &params) {
    m_bond = make_shared_from_args<CoreBond, std::string, double, int, std::vector<double>,
				   std::vector<double>, std::vector<double>>
      (params, "filename", "maxval", "noverlaps", "para_a","para_b","para_c");
  }
  
};
}
}

#endif
