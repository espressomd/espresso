/*
 * Copyright (C) 2010-2019 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef SCRIPT_INTERFACE_COLLISION_DETECTION_COLLISION_DETECTION_HPP
#define SCRIPT_INTERFACE_COLLISION_DETECTION_COLLISION_DETECTION_HPP

#include "../ScriptInterfaceBase.hpp"
#include "core/collision.hpp"

#ifdef COLLISION_DETECTION

namespace ScriptInterface {
namespace CollisionDetection {

class CollisionDetection : public AutoParameters<CollisionDetection> {
public:
  CollisionDetection() {
    add_parameters(
        {//{"bond_three_particles", collision_params.bond_three_particles},
         //{"three_particle_binding_angle_resolution",
         // collision_params.three_particle_angle_resolution},
         {"active", collision_params.active},
         {"distance", collision_params.distance},
         {"rate", collision_params.rate},
         {"particle_type", collision_params.particle_type},
         {"particle_type_after_collision", collision_params.particle_type_after_collision},
         {"vs_particle_type", collision_params.vs_particle_type},
         {"distance_vs_particle", collision_params.distance_vs_particle},
         {"bond_type", collision_params.bond_type},
         {"vs_bond_type", collision_params.vs_bond_type}});
  };
  Variant call_method(const std::string &name,
                      const VariantMap &params) override {
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
