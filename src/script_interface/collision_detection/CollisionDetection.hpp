/*
  Copyright (C) 2010,2011,2012,2013,2014 The ESPResSo project
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

#ifndef SCRIPT_INTERFACE_COLLISION_DETECTION_COLLISION_DETECTION_HPP
#define SCRIPT_INTERFACE_COLLISION_DETECTION_COLLISION_DETECTION_HPP

#include "core/collision.hpp"
#include "../ScriptInterfaceBase.hpp"

#ifdef COLLISION_DETECTION

namespace ScriptInterface {
namespace CollisionDetection {

class CollisionDetection : public AutoParameters {
public:
  CollisionDetection() {
    add_parameters({
      {"mode", collision_params.mode},
      {"exception_on_collision", collision_params.exception_on_collision},
      
      {"bond_centers",collision_params.bond_centers},
      {"bond_vs",collision_params.bond_vs},
      {"bond_three_particles",collision_params.bond_three_particles},
      {"three_particle_binding_angle_resolution",collision_params.three_particle_angle_resolution},

      {"distance",collision_params.distance},
      {"distance_glued_particle_to_vs",collision_params.dist_glued_part_to_vs},
      {"vs_placement", collision_params.vs_placement},

      {"part_type_vs",collision_params.vs_particle_type},
      {"part_type_to_be_glued",collision_params.part_type_to_be_glued},
      {"part_type_to_attach_vs_to",collision_params.part_type_to_attach_vs_to},
      {"part_type_after_glueing",collision_params.part_type_after_glueing}
    });
  };
  Variant call_method(const std::string& name, const VariantMap& params) override {
      if (name == "validate") {
          return validate_collision_parameters();
      };
    return none;
  };
};

} /* namespace CollisionDetection */
} /* namespace ScriptInterface */

#endif
#endif
