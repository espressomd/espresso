#ifndef SCRIPT_INTERFACE_WALBERLA_WALBERLABLOCKFOREST_HPP
#define SCRIPT_INTERFACE_WALBERLA_WALBERLABLOCKFOREST_HPP

#include "walberla_bridge/WalberlaBlockForest.hpp"

#include "grid.hpp"
#include "script_interface/ScriptInterface.hpp"
#include "script_interface/auto_parameters/AutoParameters.hpp"
#include "utils/Vector.hpp"

#include <cmath>

namespace ScriptInterface::walberla {

class WalberlaBlockForest : public AutoParameters<WalberlaBlockForest> {
public:
  std::shared_ptr<::walberla::WalberlaBlockForest> m_blockforest;

public:
  WalberlaBlockForest() {
    add_parameters({{"grid_dimensions", AutoParameter::read_only,
                     [this]() { return m_blockforest->get_grid_dimensions(); }},
                    {"ghost_layers", AutoParameter::read_only,
                     [this]() { return m_blockforest->get_ghost_layers(); }}});
  };

  void do_construct(VariantMap const &args) override {
    const auto agrid = get_value<double>(args, "agrid");
    const auto box_size = box_geo.length();

    if (agrid <= 0) {
      throw std::domain_error("agrid has to be >=0!");
    }

    const Utils::Vector3i grid_dimensions{
        static_cast<int>(std::round(box_size[0] / agrid)),
        static_cast<int>(std::round(box_size[1] / agrid)),
        static_cast<int>(std::round(box_size[2] / agrid))};
    for (int i : {0, 1, 2}) {
      if (std::abs(grid_dimensions[i] * agrid - box_size[i]) / box_size[i] >
          std::numeric_limits<double>::epsilon()) {
        throw std::runtime_error(
            "Box length not commensurate with agrid in direction " +
            std::to_string(i) + " length " + std::to_string(box_size[i]) +
            " agrid " + std::to_string(agrid));
      }
    }
    m_blockforest = std::make_shared<::walberla::WalberlaBlockForest>(
        grid_dimensions, node_grid, get_value<int>(args, "ghost_layers"));
  }
};

} // namespace ScriptInterface::walberla

#endif
