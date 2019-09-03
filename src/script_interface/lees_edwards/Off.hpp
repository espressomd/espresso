#ifndef SCRIPT_INTERFACE_LEES_EDWARDS_OFF_HPP
#define SCRIPT_INTERFACE_LEES_EDWARDS_OFF_HPP

#include "../ScriptInterfaceBase.hpp"
#include "config.hpp"
#include "core/lees_edwards.hpp"

#ifdef LEES_EDWARDS

namespace ScriptInterface {
namespace LeesEdwards {

class Off : public Protocol {
public:
  Off()
      : m_protocol{new ::LeesEdwards::ActiveProtocol{::LeesEdwards::Off()}} {};
  std::shared_ptr<::LeesEdwards::ActiveProtocol> protocol() override {
    return m_protocol;
  }

private:
  std::shared_ptr<::LeesEdwards::ActiveProtocol> m_protocol;
}; // Class ProtocolOff;

} // namespace LeesEdwards
} // namespace ScriptInterface

#endif
#endif
