/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
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
#ifndef VIRTUAL_SITES_VIRTUAL_SITES_HPP
#define VIRTUAL_SITES_VIRTUAL_SITES_HPP


/** \file virtual_sites.hpp
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
#include <memory>

/** @brief Base class for virtual sites implementations */
class VirtualSites {
  public:
    VirtualSites() :  m_have_velocity(true) {};
    /** @brief Update positions and/or velocities of virtual sites 

    * Velocities are only updated update_velocities() return true 
    * @param recalc_positions can be used to skip the reculation of positions 
    */
    virtual void update(bool recalc_positions=true) const =0;
    /** Back-transfer forces (and torques) to non-virtual particles */
    virtual void back_transfer_forces_and_torques() const =0;
    /** @brief Number of pressure contributions */
    virtual int n_pressure_contribs() const {return 0;};
    /** @brief Pressure contribution() */
    virtual void pressure_and_stress_tensor_contribution(double* pressure, double* stress_tensor) const {};    
    /** @brief Enable/disable velocity calculations for vs */
    void set_have_velocity(const bool& v) { m_have_velocity=v; };
    const bool& have_velocity() const { return m_have_velocity; };
    /** @brief Is a ghost communication needed after position updates */
    virtual bool need_ghost_comm_after_pos_update() const =0;
    /** Is a ghost comm needed before a velocity update */
    virtual bool need_ghost_comm_before_vel_update() const =0;
    /** Is a ghost comm needed before the back_transfer */
    virtual bool need_ghost_comm_before_back_transfer() const =0;
    private:
      bool m_have_velocity;
   };


#endif
#endif
