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
#ifndef VIRTUAL_SITES_VIRTUAL_SITES_HPP
#define VIRTUAL_SITES_VIRTUAL_SITES_HPP

/** \file
 *  This file contains routine to handle virtual sites
 *  Virtual sites are like particles, but they will be not integrated.
 *  Step performed for virtual sites:
 *  - update virtual sites
 *  - calculate forces
 *  - distribute forces
 *  - move no-virtual particles
 *  - update virtual sites
 */

#ifdef VIRTUAL_SITES
#include <utils/Vector.hpp>

/** @brief Base class for virtual sites implementations */
class VirtualSites {
public:
  VirtualSites() = default;
  virtual ~VirtualSites() = default;

  /**
   * @brief Update positions and velocities of virtual sites.
   */
  virtual void update() const {}
  /** Back-transfer forces (and torques) to non-virtual particles. */
  virtual void back_transfer_forces_and_torques() const {}
  /** @brief Called after force calculation (and before rattle/shake) */
  virtual void after_force_calc(){};
  virtual void after_lb_propagation(){};
  /** @brief Pressure contribution. */
  virtual Utils::Matrix<double, 3, 3> stress_tensor() const { return {}; };
  /** @brief Enable/disable quaternion calculations for vs.*/
  void set_have_quaternion(const bool &have_quaternion) {
    m_have_quaternion = have_quaternion;
  };
  bool get_have_quaternion() const { return m_have_quaternion; };

private:
  bool m_have_quaternion = false;
};

#endif
#endif
