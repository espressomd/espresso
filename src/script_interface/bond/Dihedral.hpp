#ifndef SCRIPT_INTERFACE_BONDS_DIHEDRAL_HPP
#define SCRIPT_INTERFACE_BONDS_DIHEDRAL_HPP

#include "Bond.hpp"

#include "core/bond/Dihedral.hpp"
#include "get_value.hpp"

#include <utility>
  
namespace ScriptInterface {
namespace Bond {
class Dihedral : public AutoParameters<Bond> {
  using CoreBond = ::Bond::Dihedral;
  std::shared_ptr<CoreBond> m_bond;

public:
 Dihedral () {
    add_parameters(
		   {{"mult",
			 [this](Variant const &v) { m_bond->mult() = get_value<double>(v); },
			 [this]() { return m_bond->mult(); }},
		       {"bend",
			   [this](Variant const &v) { m_bond->bend() = get_value<double>(v); },
			   [this]() { return m_bond->bend(); }},
			 {"phase",
			     [this](Variant const &v) { m_bond->phase() = get_value<double>(v); },
			     [this]() { return m_bond->phase(); }}
		   });
  }

  std::shared_ptr<::Bond::Bond> bond() const {
    return std::static_pointer_cast<::Bond::Bond>(m_bond);
  }

  void construct(const VariantMap &params) {
    m_bond = make_shared_from_args<CoreBond, double, double, double>(
        params, "mult", "bend", "phase");
  }
};
}
}

#endif
