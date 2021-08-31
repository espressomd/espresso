#include "config.hpp"
#include "core/lees_edwards.hpp"
#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"

namespace ScriptInterface {
namespace LeesEdwards {

class LinearShear : public Protocol {
public:
  LinearShear()
      : m_protocol{
            new ::LeesEdwards::ActiveProtocol{::LeesEdwards::LinearShear()}} {
    add_parameters(
        {{"initial_pos_offset",
          boost::get<::LeesEdwards::LinearShear>(*m_protocol)
              .m_initial_pos_offset},
         {"shear_velocity",
          boost::get<::LeesEdwards::LinearShear>(*m_protocol).m_shear_velocity},
         {"time_0",
          boost::get<::LeesEdwards::LinearShear>(*m_protocol).m_time_0}});
#ifdef LB_WALBERLA
    ::lees_edwards_active_protocol = m_protocol;
#endif
  }
  std::shared_ptr<::LeesEdwards::ActiveProtocol> protocol() override {
    return m_protocol;
  }

private:
  std::shared_ptr<::LeesEdwards::ActiveProtocol> m_protocol;
}; // Class ProtocolLinearShear;

} // namespace LeesEdwards
} // namespace ScriptInterface
