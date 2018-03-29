#ifndef SCRIPT_INTERFACE_BONDS_THERMALIZED_BOND_HPP
#define SCRIPT_INTERFACE_BONDS_THERMALIZED_BOND_HPP

#include "Bond.hpp"

#include "core/bond/ThermalizedBond.hpp"
#include "get_value.hpp"

#include <utility>
  
namespace ScriptInterface {
namespace Bond {
class ThermalizedBond : public AutoParameters<Bond> {
  using CoreBond = ::Bond::ThermalizedBond;
  std::shared_ptr<CoreBond> m_bond;

public:
  ThermalizedBond() {
    add_parameters(
		   {{"temp_com",
			 [this](Variant const &v) { m_bond->temp_com() = get_value<double>(v); },
			 [this]() { return m_bond->temp_com(); }},
		       {"gamma_com",
			   [this](Variant const &v) { m_bond->gamma_com() = get_value<double>(v); },
			   [this]() { return m_bond->gamma_com(); }},
			 {"temp_distance",
			     [this](Variant const &v) { m_bond->temp_distance() = get_value<double>
			       (v); },
			     [this]() { return m_bond->temp_distance(); }},
			   {"r_cut",
			       [this](Variant const &v) { m_bond->r_cut() = get_value<double>(v); },
			       [this]() { return m_bond->r_cut(); }},
			     {"pref1_com",
				 [this](Variant const &v) { m_bond->pref1_com() = get_value<double>
				   (v); },
				 [this]() { return m_bond->pref1_com(); }},
			       {"pref2_com",
				   [this](Variant const &v) { m_bond->pref2_com() = get_value<double>
				     (v); },
				   [this]() { return m_bond->pref2_com(); }},
				 {"pref1_dist",
				     [this](Variant const &v) { m_bond->pref1_dist() = get_value
				       <double>(v); },
				     [this]() { return m_bond->pref1_dist(); }},
				   {"pref2_dist",
				       [this](Variant const &v) { m_bond->pref2_dist() = get_value
					 <double>(v); },
				       [this]() { return m_bond->pref2_dist(); }}
		   
		   });
  }

  std::shared_ptr<::Bond::Bond> bond() const {
    return std::static_pointer_cast<::Bond::Bond>(m_bond);
  }

  void construct(const VariantMap &params) {
    m_bond = make_shared_from_args<CoreBond, double, double, double, double, double, double,
				   double, double, double>
      (params, "temp_com","gamma_com","temp_distance","gamma_distance", "r_cut", "pref1_com",
       "pref2_com", "pref1_dist", "pref2_dist");
  }
  
};
}
}

#endif
