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
#include "LBBoundary.hpp"
#include "errorhandling.hpp"
#include "grid_based_algorithms/lb_boundaries.hpp"

namespace LBBoundaries {
// Vector3d LBBoundary::get_force() {
//   double tmp_force[3];
//   int n = 0;
//   for (auto it = lbboundaries.begin(); it != lbboundaries.end(); ++it, ++n) {
//     if (&(**it) == this)
// 	break;
//   }
//   if (n == lbboundaries.size())
//     throw std::runtime_error("You probably tried to get the force of an
//     lbboundary that was not added to system.lbboundaries.");
//   lbboundary_get_force(n, tmp_force);
//   return Vector3d{tmp_force[0], tmp_force[1], tmp_force[2]};
// }
}
