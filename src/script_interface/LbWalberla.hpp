/*
  Copyright (C) 2010-2018 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
  Max-Planck-Institute for Polymer Research, Theory Group

  This file is part of ESPResSo.

  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef SCRIPT_INTERFACE_LBWALBERLA_HPP
#define SCRIPT_INTERFACE_LBWALBERLA_HPP

#include "config.hpp"

#ifdef LB_WALBERLA

#include "ScriptInterfaceBase.hpp"
#include "auto_parameters/AutoParameters.hpp"
#include "grid.hpp"

#include "core/grid_based_algorithms/LbWalberla.hpp"
#include "core/integrate.hpp"

namespace ScriptInterface {

class LbWalberla : public AutoParameters<LbWalberla> {
public:
  LbWalberla() {
    add_parameters({{"agrid", AutoParameter::read_only,
                     [this]() { return m_lb_walberla->get_grid_spacing(); }},
                    {"visc", AutoParameter::read_only,
                     [this]() { return m_lb_walberla->get_viscosity(); }}});
  }
  void construct(VariantMap const &args) override {
    m_lb_walberla = std::make_shared<::LbWalberla>(
        get_value<double>(args.at("visc")), get_value<double>(args.at("agrid")),
        Vector3d{{box_l[0], box_l[1], box_l[2]}},
        Vector3i{{node_grid[0], node_grid[1], node_grid[2]}}, skin);
  }
  Variant call_method(const std::string &name,
                      const VariantMap &params) override {
    if (name == "node_get_v") {
      auto const node = get_value<Vector3i>(params.at("node"));
      auto const v = m_lb_walberla->get_node_velocity(node);
      return v;
    }
    if (name == "node_set_v") {
      auto const node = get_value<Vector3i>(params.at("node"));
      auto const v = get_value<Vector3d>(params.at("v"));
      m_lb_walberla->set_node_velocity(node, v);
    }
    return none;
  };

private:
  std::shared_ptr<::LbWalberla> m_lb_walberla;
};

} // namespace ScriptInterface
#endif
#endif

