/*
Copyright (C) 2010-2018 The ESPResSo project

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
#ifndef VIRTUAL_SITES_VIRTUAL_SITES_OFF_HPP
#define VIRTUAL_SITES_VIRTUAL_SITES_OFF_HPP

#include "config.hpp"
#ifdef VIRTUAL_SITES
#include "VirtualSites.hpp"

/** @brief Do-nothing virtual-sites implementation */
class VirtualSitesOff : public VirtualSites {
  /** @brief Update positions and/or velocities of virtual sites

  * Velocities are only updated update_velocities() return true
  * @param recalc_positions can be used to skip the reculation of positions
  */
  void update(bool recalc_positions = true) const override{};
  /** Back-transfer forces (and torques) to non-virtual particles */
  void back_transfer_forces_and_torques() const override{};
  /** @brief Is a ghost communication needed after position updates */
  bool need_ghost_comm_after_pos_update() const override { return false; }
  /** Is a ghost comm needed before a velocity update */
  bool need_ghost_comm_before_vel_update() const override { return false; };
  /** Is a ghost comm needed before back_transfer */
  bool need_ghost_comm_before_back_transfer() const override { return false; };
};

#endif
#endif
