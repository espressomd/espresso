#ifndef SCRIPT_INTERFACE_BONDS_UMBRELLA_HPP
#define SCRIPT_INTERFACE_BONDS_UMBRELLA_HPP

#include "Bond.hpp"

#include "core/bond/Umbrella.hpp"
#include "get_value.hpp"

#include <utility>
  
namespace ScriptInterface {
namespace Bond {
class Umbrella : public AutoParameters<Bond> {
  using CoreBond = ::Bond::Umbrella;
  std::shared_ptr<CoreBond> m_bond;

public:
 Umbrella () {
    add_parameters(
        {{"k",
          [this](Variant const &v) { m_bond->k() = get_value<double>(v); },
          [this]() { return m_bond->k(); }},
{"dir",
          [this](Variant const &v) { m_bond->dir() = get_value<int>(v); },
          [this]() { return m_bond->dir(); }},
{"r",
    [this](Variant const &v) { m_bond->r() = get_value<double>(v); },
          [this]() { return m_bond->r(); }}
});
  }

  std::shared_ptr<::Bond::Bond> bond() const {
    return std::static_pointer_cast<::Bond::Bond>(m_bond);
  }

  void construct(const VariantMap &params) {
    m_bond = make_shared_from_args<CoreBond, double, int, double>(
        params, "k", "dir", "r");
  }
};
}
}

#endif
