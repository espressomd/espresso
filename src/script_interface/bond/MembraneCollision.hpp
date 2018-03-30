#ifndef SCRIPT_INTERFACE_MEMBRANE_COLLISION_HPP
#define SCRIPT_INTERFACE_MEMBRANE_COLLISION_HPP

#include "Bond.hpp"

#include "core/bond/MembraneCollision.hpp"
#include "get_value.hpp"

#include <utility>
  
namespace ScriptInterface {
namespace Bond {
class MembraneCollision : public AutoParameters<Bond> {
  using CoreBond = ::Bond::MembraneCollision;
  std::shared_ptr<CoreBond> m_bond;

public:
 MembraneCollision () {
  }

  std::shared_ptr<::Bond::Bond> bond() const {
    return std::static_pointer_cast<::Bond::Bond>(m_bond);
  }

  void construct(const VariantMap &params) {
    m_bond = make_shared_from_args<CoreBond>(params);
  }
};
}
}

#endif
