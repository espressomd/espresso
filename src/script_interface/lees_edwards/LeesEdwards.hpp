
#include "Protocol.hpp"
#include "cells.hpp"
#include "config.hpp"
#include "core/grid.hpp"
#include "core/lees_edwards.hpp"
#include "integrate.hpp"
#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"

namespace ScriptInterface {
namespace LeesEdwards {

class LeesEdwards : public AutoParameters<LeesEdwards> {
public:
  LeesEdwards() : m_protocol{nullptr} {
    add_parameters(
        {{"protocol",
          [this](Variant const &value) {
            if (is_none(value)) {
              m_protocol = nullptr;
              box_geo.lees_edwards_bc().shear_velocity = 0;
              box_geo.lees_edwards_bc().pos_offset = 0;
              box_geo.set_type(BoxType::CUBOID);
              cell_structure.set_resort_particles(Cells::RESORT_LOCAL);
              return;
            }
            m_protocol = get_value<std::shared_ptr<Protocol>>(value);
            if (m_protocol) {
              box_geo.set_type(BoxType::LEES_EDWARDS);
              ::LeesEdwards::active_protocol = m_protocol->protocol();
              ::LeesEdwards::update_pos_offset(*::LeesEdwards::active_protocol,
                                               box_geo, get_sim_time());
              ::LeesEdwards::update_shear_velocity(
                  *::LeesEdwards::active_protocol, box_geo, get_sim_time());
              cell_structure.set_resort_particles(Cells::RESORT_LOCAL);
            } else {
              throw std::runtime_error(
                  "A Lees Edwards protocol needs to be passed.");
            }
          },
          [this]() {
            if (m_protocol)
              return make_variant(m_protocol);
            return make_variant(none);
          }},
         {"shear_velocity", box_geo.lees_edwards_bc().shear_velocity},
         {"pos_offset", box_geo.lees_edwards_bc().pos_offset},
         {"shear_direction", box_geo.lees_edwards_bc().shear_direction},
         {"shear_plane_normal", box_geo.lees_edwards_bc().shear_plane_normal}});
  }

private:
  std::shared_ptr<Protocol> m_protocol;

}; // Class LeesEdwards

} // namespace LeesEdwards
} // namespace ScriptInterface
