/*
  Copyright (C) 2010-2018 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
    Max-Planck-Institute for Polymer Research, Theory Group

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
/** \file
 *
 * Boundary conditions for lattice Boltzmann fluid dynamics.
 * Header file for \ref lbboundaries.cpp.
 *
 * In the current version only simple bounce back walls are implemented. Thus
 * after the streaming step, in all wall nodes all populations are bounced
 * back from where they came from. Ulf Schiller spent a lot of time
 * working on more powerful alternatives, they are to be found in the
 * lb_testing branch of espresso until the end of 2010. Now we stripped
 * down the code to a minimum, as most of it was not sufficiently
 * understandable.
 *
 * Anyone who wants to revive these, please look into the git.
 *
 */

#ifndef LBBOUNDARIES_H
#define LBBOUNDARIES_H

#include "lbboundaries/LBBoundary.hpp"
#include "utils.hpp"

#include "utils/Span.hpp"

#include <array>

namespace LBBoundaries {
using LB_Fluid = std::array<Utils::Span<double>, 19>;

extern std::vector<std::shared_ptr<LBBoundary>> lbboundaries;
#if defined(LB_BOUNDARIES) || defined(LB_BOUNDARIES_GPU)

/** Initializes the constrains in the system.
 *  This function determines the lattice sited which belong to boundaries
 *  and marks them with a corresponding flag.
 */
void lb_init_boundaries();
void lbboundary_mindist_position(const Utils::Vector3d &pos, double *mindist,
                                 double distvec[3], int *no);

int lbboundary_get_force(int no, double *f);

void add(const std::shared_ptr<LBBoundary> &);
void remove(const std::shared_ptr<LBBoundary> &);

#endif // (LB_BOUNDARIES) || (LB_BOUNDARIES_GPU)
} // namespace LBBoundaries
#endif /* LB_BOUNDARIES_H */
