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

#ifndef VIRTUAL_SITES_CENTER_OF_MASS_HPP
#define VIRTUAL_SITES_CENTER_OF_MASS_HPP

/** \file
 *  This file contains routine to handle virtual sites at the center of mass of
 * a bunch of other particles (say, a molecule). Forces acting on this center of
 * mass are distributed back onto the constituents. The position/velocity/mass
 * of the virtual site at center of mass is calculated from the
 * positions/velocities/masses of many particles.
 *
 *  Virtual sites are like particles, but they will not be integrated.
 *  Step performed for virtual sites:
 *  - update virtual sites
 *  - calculate forces
 *  - distribute forces
 *  - move non-virtual particles
 *  - update virtual sites
 */

#include "script_interface/ObjectHandle.hpp"
#include "config/config.hpp"
#include <map>
#include <utils/Vector.hpp>
#include <utils/matrix.hpp>

/** @brief Base class for virtual sites implementations */
class VirtualSitesCenterOfMass : ObjectHandle
{
public:
  VirtualSitesCenterOfMass() = default;

  virtual ~VirtualSitesCenterOfMass() = default;

  /**
   * @brief Update positions and velocities of virtual sites.
   */
  virtual void update() const {}
  
  /** @brief Back-transfer forces to non-virtual particles. */
  virtual void back_transfer_forces() const {}
  
  /**  @brief Store (mol_id, virtual_site_particle_id) pairs */
  std::unordered_map<int, int> vitual_site_id_for_mol_id() {}

private:
  struct ComInfo {
    double total_mass = 0.0;
    Utils::Vector3d weighted_position_sum = {0., 0., 0.};
  };

  /**  @brief Store (mol_id, com_info) pairs */
  std::unordered_map<int, std::shared_ptr<ComInfo>> com_by_mol_id;

  bool m_have_quaternion = false;
  bool m_override_cutoff_check = false;
};

#endif
