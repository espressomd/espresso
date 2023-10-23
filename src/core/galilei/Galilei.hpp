/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#pragma once

#include "system/System.hpp"

#include <utils/Vector.hpp>

class Galilei {
public:
  /** Stop particle motion by setting the velocity of each particle to zero.
   *  @param system the system to modify
   *  @param omega  if true, also set particle angular velocities to zero
   */
  void kill_particle_motion(System::System &system, bool omega) const;

  /** Set all the forces acting on the particles to zero.
   *  @param system  the system to modify
   *  @param torque  if true, also set particle torques to zero
   */
  void kill_particle_forces(System::System &system, bool torque) const;

  /** Calculate the CMS of the system */
  Utils::Vector3d calc_system_cms_position(System::System const &system) const;

  /** Calculate the CMS velocity of the system */
  Utils::Vector3d calc_system_cms_velocity(System::System const &system) const;

  /** Remove the CMS velocity */
  void galilei_transform(System::System &system) const;
};
