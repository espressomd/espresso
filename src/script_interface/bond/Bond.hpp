#ifndef SCRIPT_INTERFACE_BOND_BOND_HPP
#define SCRIPT_INTERFACE_BOND_BOND_HPP

#include "ScriptInterfaceBase.hpp"

#include "core/bond/Bond.hpp"

#include <memory>

namespace ScriptInterface {
  namespace Bond {
    class Bond : public ScriptInterfaceBase {
      int m_id = -1;
    public:
      using CoreBond = ::Bond::Bond;

      int &id() { return m_id; }
      
      virtual std::shared_ptr<CoreBond> bond() const = 0;
    };
  }}

#endif
