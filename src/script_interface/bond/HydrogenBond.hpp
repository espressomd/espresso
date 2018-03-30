#ifndef SCRIPT_INTERFACE_BONDS_HYDROGEN_BOND_HPP
#define SCRIPT_INTERFACE_BONDS_HYDROGEN_BOND_HPP

#include "Bond.hpp"

#include "core/bond/HydrogenBond.hpp"
#include "get_value.hpp"

#include <utility>
  
namespace ScriptInterface {
namespace Bond {
class HydrogenBond : public AutoParameters<Bond> {
  using CoreBond = ::Bond::HydrogenBond;
  std::shared_ptr<CoreBond> m_bond;

public:
  HydrogenBond() {
    add_parameters(
		   {{"r0",
			 [this](Variant const &v) { m_bond->r0() = get_value<double>(v); },
			 [this]() { return m_bond->r0(); }},
		       {"alpha",
			   [this](Variant const &v) { m_bond->alpha() = get_value<double>(v); },
			   [this]() { return m_bond->alpha(); }},
			 {"E0",
			     [this](Variant const &v) { m_bond->E0() = get_value<double>(v); },
			     [this]() { return m_bond->E0(); }},
			   {"kd",
			       [this](Variant const &v) { m_bond->kd() = get_value<double>(v); },
			       [this]() { return m_bond->kd(); }},
			     {"sigma1",
				 [this](Variant const &v) { m_bond->sigma1() = get_value<double>
				   (v); },
				 [this]() { return m_bond->sigma1(); }},
			       {"sigma2",
				   [this](Variant const &v) { m_bond->sigma2() = get_value<double>
				     (v); },
				   [this]() { return m_bond->sigma2(); }},
				 {"psi10",
				     [this](Variant const &v) { m_bond->psi10() = get_value<double>
				       (v); },
				     [this]() { return m_bond->psi10(); }},
				   {"psi20",
				       [this](Variant const &v) { m_bond->psi20() = get_value
					 <double>
					 (v); },
				       [this]() { return m_bond->psi20(); }},
				     {"E0sb",
					 [this](Variant const &v) { m_bond->E0sb() = get_value
					   <double>
					   (v); },
					 [this]() { return m_bond->E0sb(); }},
				       {"r0sb",
					   [this](Variant const &v) { m_bond->r0sb() = get_value
					     <double>
					     (v); },
					   [this]() { return m_bond->r0sb(); }},
				       {"alphasb",
					   [this](Variant const &v) { m_bond->alphasb() =
					     get_value<double>
					     (v); },
					   [this]() { return m_bond->alphasb(); }},
					 {"f2",
					     [this](Variant const &v) { m_bond->f2() = get_value
					       <double>
					       (v); },
					     [this]() { return m_bond->f2(); }},
					   {"f3",
					       [this](Variant const &v) { m_bond->f3() = 
						 get_value <double>
						 (v); },
					       [this]() { return m_bond->f3(); }},
		   });
  }

  std::shared_ptr<::Bond::Bond> bond() const {
    return std::static_pointer_cast<::Bond::Bond>(m_bond);
  }

  void construct(const VariantMap &params) {
    m_bond = make_shared_from_args<CoreBond, double, double, double, double, double,double, double,
				   double, double, double, double, double, double>
      (params, "r0","alpha", "E0", "kd", "sigma1", "sigma2", "psi10", "psi20", "E0sb", "r0sb",
       "alphasb", "f2", "f3");
  }
  
};
}
}

#endif
