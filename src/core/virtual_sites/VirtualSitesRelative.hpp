/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
  Copyright (C) 2010,2011 Rudolf Weeber

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
#ifndef VIRTUAL_SITES_VIRTUAL_SITES_RELATIVE_HPP
#define VIRTUAL_SITES_VIRTUAL_SITES_RELATIVE_HPP

#include "config.hpp" 
#ifdef VIRTUAL_SITES_RELATIVE

#include "virtual_sites.hpp"
#include "particle_data.hpp"
#include "communication.hpp"


   /** @brief Virtual sites implementation for rigid bodies */
   class VirtualSitesRelative : public VirtualSites {
    public:
    VirtualSitesRelative() {};
    /** @brief Update positions and/or velocities of virtual sites 

    * Velocities are only updated have_velocity() return true 
    * @param recalc_positions can be used to skip the reculation of positions 
    */
    void update(bool recalc_positions=true) const override;
    /** Back-transfer forces (and torques) to non-virtual particles */
    void back_transfer_forces_and_torques() const override;
    /** @brief Is a ghost communication needed before position updates */
    bool need_ghost_comm_after_pos_update() const override { return true;} 
    /** Is a ghost comm needed before a velocity update */
    bool need_ghost_comm_before_vel_update() const override {return (n_nodes>1) && have_velocity();};
    bool need_ghost_comm_before_back_transfer() const override {return true;};
    int n_pressure_contribs() const override {return 1;};
    void pressure_and_stress_tensor_contribution(double* pressure, double* stress_tensor) const override;    
    
    private:
    void update_pos(Particle& p) const;
    void update_vel(Particle& p) const;
   };

#endif

#endif
