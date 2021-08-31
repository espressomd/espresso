#ifndef SCRIPT_INTERFACE_LEES_EDWARDS_OFF_HPP
#define SCRIPT_INTERFACE_LEES_EDWARDS_OFF_HPP

#include "config.hpp"
#include "core/lees_edwards.hpp"
#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"

namespace ScriptInterface {
namespace LeesEdwards {

class Off : public Protocol {
public:
  Off() : m_protocol{new ::LeesEdwards::ActiveProtocol{::LeesEdwards::Off()}} {
#ifdef LB_WALBERLA
    ::lees_edwards_active_protocol = m_protocol;
#endif
  }
  std::shared_ptr<::LeesEdwards::ActiveProtocol> protocol() override {
    return m_protocol;
  }

private:
  std::shared_ptr<::LeesEdwards::ActiveProtocol> m_protocol;
}; // Class ProtocolOff;

} // namespace LeesEdwards
} // namespace ScriptInterface

#endif
