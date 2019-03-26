/*
  Copyright (C) 2010-2018 The ESPResSo project
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

#include "communication.hpp"
#include "particle_data.hpp"
#include "virtual_sites.hpp"

/** @brief Virtual sites implementation for rigid bodies */
class VirtualSitesRelative : public VirtualSites {
public:
  VirtualSitesRelative() = default;
  /** @copydoc VirtualSites::update */
  void update(bool recalc_positions) const override;
  /** @copydoc VirtualSites::back_transfer_forces_and_torques */
  void back_transfer_forces_and_torques() const override;
  /** @copydoc VirtualSites::need_ghost_comm_after_pos_update */
  bool need_ghost_comm_after_pos_update() const override { return true; }
  /** @copydoc VirtualSites::need_ghost_comm_before_vel_update */
  bool need_ghost_comm_before_vel_update() const override {
    return (n_nodes > 1) && get_have_velocity();
  };
  /** @copydoc VirtualSites::need_ghost_comm_before_back_transfer */
  bool need_ghost_comm_before_back_transfer() const override { return true; };
  /** @copydoc VirtualSites::n_pressure_contribs */
  int n_pressure_contribs() const override { return 1; };
  /** @copydoc VirtualSites::pressure_and_stress_tensor_contribution */
  void
  pressure_and_stress_tensor_contribution(double *pressure,
                                          double *stress_tensor) const override;

private:
  void update_pos(Particle &p) const;
  void update_vel(Particle &p) const;
  /** @brief Update the orientation of the virtual particles with respect to the
   * real particle.
   */
  void update_virtual_particle_quaternion(Particle &p) const;
};

#endif

#endif
