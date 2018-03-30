#ifndef SCRIPT_INTERFACE_BONDS_ANGLE_DIST_HPP
#define SCRIPT_INTERFACE_BONDS_ANGLE_DIST_HPP

#include "Bond.hpp"

#include "core/bond/AngleDist.hpp"
#include "get_value.hpp"

#include <utility>
  
namespace ScriptInterface {
namespace Bond {
class AngleDist : public AutoParameters<Bond> {
  using CoreBond = ::Bond::AngleDist;
  std::shared_ptr<CoreBond> m_bond;

public:
  AngleDist() {
    add_parameters(
		   {{"bend",
			 [this](Variant const &v) { m_bond->bend() = get_value<double>(v); },
			 [this]() { return m_bond->bend(); }},
		       {"phimin",
			   [this](Variant const &v) { m_bond->phimin() = get_value<double>(v); },
			   [this]() { return m_bond->phimin(); }},
			 {"distmin",
			     [this](Variant const &v) { m_bond->distmin() = get_value<double>(v); },
			     [this]() { return m_bond->distmin(); }},
			   {"phimax",
			       [this](Variant const &v) { m_bond->phimax() = get_value<double>(v); },
			       [this]() { return m_bond->phimax(); }},
			     {"distmax",
				 [this](Variant const &v) { m_bond->distmax() = get_value<double>
				   (v); },
				 [this]() { return m_bond->distmax(); }},
		   });
  }

  std::shared_ptr<::Bond::Bond> bond() const {
    return std::static_pointer_cast<::Bond::Bond>(m_bond);
  }

  void construct(const VariantMap &params) {
    m_bond = make_shared_from_args<CoreBond, double, double, double, double, double>
      (params, "bend", "phimin", "distmin", "phimax", "distmax");
  }
  
};
}
}

#endif
