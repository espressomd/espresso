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
#ifndef VIRTUAL_SITES_VIRTUAL_SITES_INERTIALESS_TRACERS_HPP
#define VIRTUAL_SITES_VIRTUAL_SITES_INERTIALESS_TRACERS_HPP

#include "config.hpp"
#ifdef VIRTUAL_SITES
#include "VirtualSites.hpp"

#ifdef VIRTUAL_SITES_INERTIALESS_TRACERS
/** @brief Virtual sites which are advected with an lb fluid. Forces on them are
 * instantaneously transferred to the fluid
 */
class VirtualSitesInertialessTracers : public VirtualSites {
  void update(bool recalc_positions) const override{};
  void back_transfer_forces_and_torques() const override{};
  void after_force_calc() override;
  void after_lb_propagation() override;
  bool need_ghost_comm_after_pos_update() const override { return false; }
  bool need_ghost_comm_before_vel_update() const override { return false; };
  bool need_ghost_comm_before_back_transfer() const override { return false; };
};

#endif
#endif
#endif
