/*
 * Copyright (C) 2010-2022 The ESPResSo project
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
#ifndef VIRTUAL_SITES_VIRTUAL_SITES_INERTIALESS_TRACERS_HPP
#define VIRTUAL_SITES_VIRTUAL_SITES_INERTIALESS_TRACERS_HPP

#include "config/config.hpp"

#ifdef VIRTUAL_SITES_INERTIALESS_TRACERS

#include "VirtualSites.hpp"

/** @brief Virtual sites which are advected with an lb fluid. Forces on them are
 * instantaneously transferred to the fluid
 */
class VirtualSitesInertialessTracers : public VirtualSites {
  void after_force_calc(double time_step) override;
  void after_lb_propagation(double time_step) override;
};

#endif // VIRTUAL_SITES_INERTIALESS_TRACERS
#endif
