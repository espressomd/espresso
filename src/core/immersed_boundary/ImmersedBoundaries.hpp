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

#pragma once

#include "config/config.hpp"

#include "cell_system/CellStructure.hpp"
#include "system/Leaf.hpp"

#include <cassert>
#include <cstddef>
#include <vector>

struct IBMVolCons;

class ImmersedBoundaries : public System::Leaf<ImmersedBoundaries> {
public:
  ImmersedBoundaries() : VolumeInitDone(false), BoundariesFound(false) {
    VolumesCurrent.resize(10);
  }
  void init_volume_conservation(CellStructure &cs);
  void volume_conservation(CellStructure &cs);
  void register_softID(IBMVolCons &bond);

private:
  void calc_volumes(CellStructure &cs);
  void calc_volume_force(CellStructure &cs);

  std::vector<double> VolumesCurrent;
  bool VolumeInitDone;
  bool BoundariesFound;
};
