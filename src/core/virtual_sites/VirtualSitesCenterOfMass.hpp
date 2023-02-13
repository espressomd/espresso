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
 *  This file contains routine to handle virtual sites at the center of mass of a bunch of other particles (say, a molecule).
 *  Forces acting on this center of mass are distributed back onto the constituents.
 *  The position/velocity/mass of the virtual site at center of mass is calculated from the positions/velocities/masses of many particles.
 * 
 *  Virtual sites are like particles, but they will not be integrated.
 *  Step performed for virtual sites:
 *  - update virtual sites
 *  - calculate forces
 *  - distribute forces
 *  - move non-virtual particles
 *  - update virtual sites
 */

#include <map>
#include "config/config.hpp"
#include <utils/Vector.hpp>
#include <utils/matrix.hpp>

/** @brief Base class for virtual sites implementations */
class VirtualSitesCenterOfMass {
public:
  VirtualSitesCenterOfMass() = default;
  virtual ~VirtualSitesCenterOfMass() = default;

  /**
   * @brief Update positions and velocities of virtual sites.
   */
  virtual void update() const {}
  /** Back-transfer forces (and torques) to non-virtual particles. */
  virtual void back_transfer_forces_and_torques() const {}
  /** @brief Called after force calculation (and before rattle/shake) */
  virtual void after_force_calc() {}
  virtual void after_lb_propagation(double) {}
  /** @brief Pressure contribution. */
  virtual Utils::Matrix<double, 3, 3> pressure_tensor() const { return {}; }
  /** @brief Enable/disable quaternion calculations for vs.*/
  void set_have_quaternion(const bool &have_quaternion) {
    m_have_quaternion = have_quaternion;
  }
  bool have_quaternions() const { return m_have_quaternion; }
  /**  @brief Enable/disable override for the vs cutoff check */
  void set_override_cutoff_check(const bool &override_cutoff_check) {
    m_override_cutoff_check = override_cutoff_check;
  }
  bool get_override_cutoff_check() const { return m_override_cutoff_check; }
  /**  @brief Store (mol_id, virtual_site_particle_id) pairs */
  std::map<int,int> vitual_site_id_for_mol_id() {}  

private:
  bool m_have_quaternion = false;
  bool m_override_cutoff_check = false;
};

#endif
