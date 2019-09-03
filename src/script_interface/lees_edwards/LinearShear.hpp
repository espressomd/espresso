#include "../ScriptInterfaceBase.hpp"
#include "config.hpp"
#include "core/lees_edwards.hpp"

#ifdef LEES_EDWARDS

namespace ScriptInterface {
namespace LeesEdwards {

class LinearShear : public Protocol {
public:
  LinearShear()
      : m_protocol{
            new ::LeesEdwards::ActiveProtocol{::LeesEdwards::LinearShear()}} {
    add_parameters(
        {{"shear_direction", boost::get<::LeesEdwards::LinearShear>(*m_protocol)
                                 .m_shear_dir_coord},
         {"shear_plane_normal",
          boost::get<::LeesEdwards::LinearShear>(*m_protocol)
              .m_shear_plane_normal_coord},
         {"initial_pos_offset",
          boost::get<::LeesEdwards::LinearShear>(*m_protocol)
              .m_initial_pos_offset},
         {"shear_velocity", boost::get<::LeesEdwards::LinearShear>(*m_protocol)
                                .m_shear_velocity}});
  }
  std::shared_ptr<::LeesEdwards::ActiveProtocol> protocol() override {
    return m_protocol;
  }

private:
  std::shared_ptr<::LeesEdwards::ActiveProtocol> m_protocol;
}; // Class ProtocolLinearShear;

} // namespace LeesEdwards
} // namespace ScriptInterface

#endif
