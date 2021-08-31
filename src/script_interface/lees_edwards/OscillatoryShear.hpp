#include "config.hpp"
#include "core/lees_edwards.hpp"
#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"

namespace ScriptInterface {
namespace LeesEdwards {

class OscillatoryShear : public Protocol {
public:
  OscillatoryShear()
      : m_protocol{new ::LeesEdwards::ActiveProtocol{
            ::LeesEdwards::OscillatoryShear()}} {
    add_parameters(
        {{"amplitude",
          boost::get<::LeesEdwards::OscillatoryShear>(*m_protocol).m_amplitude},
         {"omega",
          boost::get<::LeesEdwards::OscillatoryShear>(*m_protocol).m_omega},
         {"time_0",
          boost::get<::LeesEdwards::OscillatoryShear>(*m_protocol).m_time_0}});
#ifdef LB_WALBERLA
    ::lees_edwards_active_protocol = m_protocol;
#endif
  }
  std::shared_ptr<::LeesEdwards::ActiveProtocol> protocol() override {
    return m_protocol;
  }

private:
  std::shared_ptr<::LeesEdwards::ActiveProtocol> m_protocol;
}; // Class ProtocolOscillatoryShear;

} // namespace LeesEdwards
} // namespace ScriptInterface
