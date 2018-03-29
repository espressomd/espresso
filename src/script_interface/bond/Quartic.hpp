#ifndef SCRIPT_INTERFACE_BONDS_QUARTIC_HPP
#define SCRIPT_INTERFACE_BONDS_QUARTIC_HPP

#include "Bond.hpp"

#include "core/bond/Quartic.hpp"
#include "get_value.hpp"

#include <utility>
  
namespace ScriptInterface {
namespace Bond {
class Quartic : public AutoParameters<Bond> {
  using CoreBond = ::Bond::Quartic;
  std::shared_ptr<CoreBond> m_bond;

public:
  Quartic() {
    add_parameters(
        {{"k0",
          [this](Variant const &v) { m_bond->k0() = get_value<double>(v); },
          [this]() { return m_bond->k0(); }},
	    {"k1",
		[this](Variant const &v) { m_bond->k1() = get_value<double>(v); },
		[this]() { return m_bond->k1(); }},
	      {"r",
		  [this](Variant const &v) { m_bond->r() = get_value<double>(v); },
		  [this]() { return m_bond->r(); }},
		{"r_cut",
		  [this](Variant const &v) { m_bond->r_cut() = get_value<double>(v); },
		  [this]() { return m_bond->r_cut(); }},
}
		   
);
  }

  std::shared_ptr<::Bond::Bond> bond() const {
    return std::static_pointer_cast<::Bond::Bond>(m_bond);
  }

  void construct(const VariantMap &params) {
    m_bond = make_shared_from_args<CoreBond, double , double, double ,double>(params, "k0",
									      "k1", "r", "r_cut");
  }
};
}
}

#endif
