#ifndef SCRIPT_INTERFACE_OIF_LOCAL_FORCES_HPP
#define SCRIPT_INTERFACE_OIF_LOCAL_FORCES_HPP

#include "Bond.hpp"

#include "core/bond/OifLocalForces.hpp"
#include "get_value.hpp"

#include <utility>
  
namespace ScriptInterface {
namespace Bond {
class OifLocalForces : public AutoParameters<Bond> {
  using CoreBond = ::Bond::OifLocalForces;
  std::shared_ptr<CoreBond> m_bond;

public:
  OifLocalForces () {
    add_parameters(
		   {{"phi0",
			 [this](Variant const &v) { m_bond->phi0() = get_value<double>(v); },
			 [this]() { return m_bond->phi0(); }},
		       {"kb",
			   [this](Variant const &v) { m_bond->kb() = get_value<double>(v); },
			   [this]() { return m_bond->kb(); }},
			 {"r0",
			     [this](Variant const &v) { m_bond->r0() = get_value<double>(v); },
			     [this]() { return m_bond->r0(); }},
			   {"ks",
			       [this](Variant const &v) { m_bond->ks() = get_value<double>(v); },
			       [this]() { return m_bond->ks(); }},
			     {"kslin",
				 [this](Variant const &v) { m_bond->kslin() = get_value<double>(v);},
				 [this]() { return m_bond->kslin(); }},
			       {"A01",
				   [this](Variant const &v) { m_bond->A01() = get_value<double>
				     (v);},
				   [this]() { return m_bond->A01(); }},
				 {"A02",
				     [this](Variant const &v) { m_bond->A02() = get_value<double>
				       (v); },
				     [this]() { return m_bond->A02(); }},
				   {"kal",
				       [this](Variant const &v) { m_bond->kal() = get_value<double>
					 (v); },
				       [this]() { return m_bond->kal(); }},
				     });
  }

  std::shared_ptr<::Bond::Bond> bond() const {
    return std::static_pointer_cast<::Bond::Bond>(m_bond);
  }

  void construct(const VariantMap &params) {
    m_bond = make_shared_from_args<CoreBond, double, double, double, double,double, double,
				   double, double>(params, "phi0", "kb", "r0", "ks",
						   "kslin","A01", "A02","kal");
  }
};
}
}

#endif
