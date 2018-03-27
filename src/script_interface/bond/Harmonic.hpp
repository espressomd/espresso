#ifndef SCRIPT_INTERFACE_BONDS_HARMONIC_HPP
#define SCRIPT_INTERFACE_BONDS_HARMONIC_HPP

#include "Bond.hpp"

#include "core/bond/Harmonic.hpp"
#include "get_value.hpp"

#include <utility>
  
namespace ScriptInterface {
namespace Bond {
class Harmonic : public AutoParameters<Bond> {
  using CoreBond = ::Bond::Harmonic;
  std::shared_ptr<CoreBond> m_bond;

public:
 Harmonic () {
    add_parameters(
        {{"k",
          [this](Variant const &v) { m_bond->k() = get_value<double>(v); },
          [this]() { return m_bond->k(); }},
{"r_0",
          [this](Variant const &v) { m_bond->r() = get_value<double>(v); },
          [this]() { return m_bond->r(); }},
{"r_cut",
    [this](Variant const &v) { m_bond->r_cut() = get_value<double>(v); },
          [this]() { return m_bond->r_cut(); }}
});
  }

  std::shared_ptr<::Bond::Bond> bond() const {
    return std::static_pointer_cast<::Bond::Bond>(m_bond);
  }

  void construct(const VariantMap &params) {
    m_bond = make_shared_from_args<CoreBond, double, double, double>(
        params, "k", "r_0", "r_cut");
  }
};
}
}

#endif
