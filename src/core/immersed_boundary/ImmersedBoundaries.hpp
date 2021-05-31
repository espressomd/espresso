/*
 * Copyright (C) 2010-2019 The ESPResSo project
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
#ifndef IMMERSED_BOUNDARY_IMMERSED_BOUNDARIES_HPP
#define IMMERSED_BOUNDARY_IMMERSED_BOUNDARIES_HPP

#include "config.hpp"

#include "CellStructure.hpp"

#include <vector>

class ImmersedBoundaries {
public:
  ImmersedBoundaries() : VolumeInitDone(false), BoundariesFound(false) {
    VolumesCurrent.resize(IBM_MAX_NUM);
  }
  void init_volume_conservation(CellStructure &cs);
  void volume_conservation(CellStructure &cs);

private:
  void calc_volumes(CellStructure &cs);
  void calc_volume_force(CellStructure &cs);

  std::vector<double> VolumesCurrent;
  bool VolumeInitDone;
  bool BoundariesFound;
};

#endif
