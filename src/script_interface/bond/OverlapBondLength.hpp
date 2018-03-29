#ifndef SCRIPT_INTERFACE_BONDS_OVERLAP_BOND_LENGTH_HPP
#define SCRIPT_INTERFACE_BONDS_OVERLAP_BOND_LENGTH_HPP

#include "Bond.hpp"

#include "core/bond/OverlapBondLength.hpp"
#include "get_value.hpp"

#include <utility>
  
namespace ScriptInterface {
namespace Bond {
class OverlapBondLength : public AutoParameters<Bond> {
  using CoreBond = ::Bond::OverlapBondLength;
  std::shared_ptr<CoreBond> m_bond;

public:
  OverlapBondLength() {
    add_parameters(
		   {{"filename",
			 [this](Variant const &v) {*m_bond->m_filename = get_value<char>(v); },
			 [this]() { return *m_bond->m_filename; }},
		       {"maxval",
			   [this](Variant const &v) { m_bond->maxval() = get_value<double>(v); },
			   [this]() { return m_bond->maxval(); }},		       
			 {"noverlaps",
			     [this](Variant const &v) { m_bond->noverlaps() = get_value<int>(v); },
			     [this]() { return m_bond->noverlaps(); }},
			   {"para_a",
			       [this](Variant const &v) { *m_bond->m_para_a = get_value<double>(v);},
			       [this]() { return m_bond->para_a(); }},
			     {"para_b",
				 [this](Variant const &v) { *m_bond->m_para_b = get_value<double>
				   (v);},
				 [this]() { return m_bond->para_b(); }},
			       {"para_c",
				   [this](Variant const &v) { *m_bond->m_para_c = get_value<double>
				     (v); },
				   [this]() { return m_bond->para_c(); }}
		      
		   });
  }

  std::shared_ptr<::Bond::Bond> bond() const {
    return std::static_pointer_cast<::Bond::Bond>(m_bond);
  }

  void construct(const VariantMap &params) {
    m_bond = make_shared_from_args<CoreBond, char, double, int, double, double, double>
      (params, "filename", "maxval", "noverlaps", "para_a","para_b","para_c");
  }
  
};
}
}

#endif
