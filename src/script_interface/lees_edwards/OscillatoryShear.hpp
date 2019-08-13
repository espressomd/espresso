#include "../ScriptInterfaceBase.hpp"
#include "config.hpp"
#include "core/lees_edwards.hpp"

#ifdef LEES_EDWARDS

namespace ScriptInterface {
namespace LeesEdwards {

class OscillatoryShear : public Protocol {
public:
  OscillatoryShear()
      : m_protocol{new ::LeesEdwards::ActiveProtocol{
            ::LeesEdwards::OscillatoryShear()}} {
    add_parameters(
        {{"shear_direction",
          boost::get<::LeesEdwards::OscillatoryShear>(*m_protocol)
              .m_shear_dir_coord},
         {"shear_plane_normal",
          boost::get<::LeesEdwards::OscillatoryShear>(*m_protocol)
              .m_shear_plane_normal_coord},
         {"amplitude",
          boost::get<::LeesEdwards::OscillatoryShear>(*m_protocol).m_amplitude},
         {"frequency",
          boost::get<::LeesEdwards::OscillatoryShear>(*m_protocol).m_frequency},
         {"time_0",
          boost::get<::LeesEdwards::OscillatoryShear>(*m_protocol).m_time_0}});
  }
  std::shared_ptr<::LeesEdwards::ActiveProtocol> protocol() override {
    return m_protocol;
  }

private:
  std::shared_ptr<::LeesEdwards::ActiveProtocol> m_protocol;
}; // Class ProtocolOscillatoryShear;

} // namespace LeesEdwards
} // namespace ScriptInterface

#endif
