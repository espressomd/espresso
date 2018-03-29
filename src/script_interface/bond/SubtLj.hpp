#ifndef SCRIPT_INTERFACE_BONDS_SUBT_LJ_HPP
#define SCRIPT_INTERFACE_BONDS_SUBT_LJ_HPP

#include "Bond.hpp"

#include "core/bond/SubtLj.hpp"
#include "get_value.hpp"

#include <utility>
  
namespace ScriptInterface {
namespace Bond {
class SubtLj : public AutoParameters<Bond> {
  using CoreBond = ::Bond::SubtLj;
  std::shared_ptr<CoreBond> m_bond;

public:
  SubtLj() {
    //params are in ia params
    //-> where do we have to set them here?
    //=> we only need to know the bond id for params
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
