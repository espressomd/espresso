#ifndef SCRIPT_INTERFACE_BONDS_IBM_TRIEL_HPP
#define SCRIPT_INTERFACE_BONDS_IBM_TRIEL_HPP

#include "Bond.hpp"

#include "core/bond/IbmTriel.hpp"
#include "get_value.hpp"

#include <utility>
  
namespace ScriptInterface {
namespace Bond {
class IbmTriel : public AutoParameters<Bond> {
  using CoreBond = ::Bond::IbmTriel;
  std::shared_ptr<CoreBond> m_bond;

public:
  IbmTriel() {
    add_parameters(
		   {{"a1",
			 [this](Variant const &v) { m_bond->a1() = get_value<double>(v); },
			 [this]() { return m_bond->a1(); }},
		       {"a2",
			   [this](Variant const &v) { m_bond->a2() = get_value<double>(v); },
			   [this]() { return m_bond->a2(); }},
			 {"b1",
			     [this](Variant const &v) { m_bond->b1() = get_value<double>(v); },
			     [this]() { return m_bond->b1(); }},
			   {"b2",
			       [this](Variant const &v) { m_bond->b2() = get_value<double>(v); },
			       [this]() { return m_bond->b2(); }},
			     {"l0",
				 [this](Variant const &v) { m_bond->l0() = get_value<double>
				   (v); },
				 [this]() { return m_bond->l0(); }},
			       {"lp0",
				   [this](Variant const &v) { m_bond->lp0() = get_value<double>
				     (v); },
				   [this]() { return m_bond->lp0(); }},
				 {"sinPhi0",
				     [this](Variant const &v) { m_bond->sinPhi0() = get_value<double>
				       (v); },
				     [this]() { return m_bond->sinPhi0(); }},
				   {"cosPhi0",
				       [this](Variant const &v) { m_bond->cosPhi0() = get_value
					 <double>
					 (v); },
				       [this]() { return m_bond->cosPhi0(); }},
				     {"area0",
					 [this](Variant const &v) { m_bond->area0() = get_value
					   <double>
					   (v); },
					 [this]() { return m_bond->area0(); }},
				       {"maxdist",
					   [this](Variant const &v) { m_bond->maxdist() = get_value
					     <double>
					     (v); },
					   [this]() { return m_bond->maxdist(); }},
				       {"elasticLaw",
					   [this](Variant const &v) { m_bond->elasticLaw() =
					     get_value<double>
					     (v); },
					   [this]() { return m_bond->elasticLaw(); }},
					 {"k1",
					     [this](Variant const &v) { m_bond->k1() = get_value
					       <double>
					       (v); },
					     [this]() { return m_bond->k1(); }},
					   {"k2",
					       [this](Variant const &v) { m_bond->k2() = 
						 get_value <double>
						 (v); },
					       [this]() { return m_bond->k2(); }},
		   });
  }

  std::shared_ptr<::Bond::Bond> bond() const {
    return std::static_pointer_cast<::Bond::Bond>(m_bond);
  }

  void construct(const VariantMap &params) {
    m_bond = make_shared_from_args<CoreBond, double, double, double, double, double,double, double,
				   double, double, double, double, double, double>
      (params, "a1","a2", "b1", "b2", "l0", "lp0", "sinPhi0", "cosPhi0", "area0", "maxdist",
       "elasticLaw", "k1", "k2");
  }
  
};
}
}

#endif
