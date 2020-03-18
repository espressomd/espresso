/*
 * Copyright (C) 2010-2019 The ESPResSo project
 * Copyright (C) 2010,2011 Rudolf Weeber
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
#ifndef VIRTUAL_SITES_VIRTUAL_SITES_RELATIVE_HPP
#define VIRTUAL_SITES_VIRTUAL_SITES_RELATIVE_HPP

#include "config.hpp"
#ifdef VIRTUAL_SITES_RELATIVE

#include "VirtualSites.hpp"

#include <utils/Vector.hpp>

/** @brief Virtual sites implementation for rigid bodies */
class VirtualSitesRelative : public VirtualSites {
public:
  VirtualSitesRelative() = default;
  /** @copydoc VirtualSites::update */
  void update() const override;
  /** @copydoc VirtualSites::back_transfer_forces_and_torques */
  void back_transfer_forces_and_torques() const override;
  /** @copydoc VirtualSites::stress_tensor */
  Utils::Matrix<double, 3, 3> stress_tensor() const override;
};

#endif

#endif
