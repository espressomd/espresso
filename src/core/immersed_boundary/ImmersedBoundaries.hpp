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
#ifndef IMMERSED_BOUNDARY_IMMERSED_BOUNDARIES_HPP
#define IMMERSED_BOUNDARY_IMMERSED_BOUNDARIES_HPP

#include "config.hpp"
#include <vector>

#ifdef IMMERSED_BOUNDARY
class ImmersedBoundaries {
public:
  ImmersedBoundaries() : MaxNumIBM(1000), VolumeInitDone(false) {
    VolumesCurrent.resize(MaxNumIBM);
  }
  void init_volume_conservation();
  void volume_conservation();
  int volume_conservation_reset_params(int bond_type, double volRef);
  int volume_conservation_set_params(int bond_type, int softID, double kappaV);
  void calc_volumes();
  void calc_volume_force();

private:
  const int MaxNumIBM;
  std::vector<double> VolumesCurrent;
  bool VolumeInitDone = false;
  bool BoundariesFound = false;
};

#endif
#endif
