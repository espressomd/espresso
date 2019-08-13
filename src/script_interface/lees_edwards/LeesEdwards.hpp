
#include "../ScriptInterfaceBase.hpp"
#include "Protocol.hpp"
#include "config.hpp"
#include "core/grid.hpp"
#include "core/lees_edwards.hpp"

#ifdef LEES_EDWARDS

namespace ScriptInterface {
namespace LeesEdwards {

class LeesEdwards : public AutoParameters<LeesEdwards> {
public:
  LeesEdwards() : m_protocol{nullptr} {
    add_parameters(
        {{"protocol",
          [this](Variant const &value) {
            m_protocol = get_value<std::shared_ptr<Protocol>>(value);
            if (m_protocol) {
              box_geo.lees_edwards_protocol = m_protocol->protocol();
              box_geo.lees_edwards_state.update(box_geo.lees_edwards_protocol,
                                                sim_time, sim_time);
            } else {
              throw std::runtime_error(
                  "A Lees Edwards protocol needs to be passed.");
            }
          },
          [this]() {
            return (m_protocol != nullptr) ? m_protocol->id() : ObjectId();
          }},
         {"shear_velocity", box_geo.lees_edwards_state.shear_velocity},
         {"pos_offset", box_geo.lees_edwards_state.pos_offset}});
  }

private:
  std::shared_ptr<Protocol> m_protocol;

}; // Class LeesEdwards

} // namespace LeesEdwards
} // namespace ScriptInterface

#endif
