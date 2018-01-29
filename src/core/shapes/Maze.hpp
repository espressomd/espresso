/*
  Copyright (C) 2010,2011,2012,2013,2014 The ESPResSo project
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

#ifndef __SHAPES_MAZE_HPP
#define __SHAPES_MAZE_HPP

#include "Shape.hpp"

namespace Shapes {
class Maze : public Shape {
public:
  Maze() : m_nsphere(0), m_dim(0), m_sphrad(0), m_cylrad(0) {}

  int calculate_dist(const double *ppos, double *dist,
                     double *vec) const override;

  int &nsphere() { return m_nsphere; }
  double &dim() { return m_dim; }
  double &sphrad() { return m_sphrad; }
  double &cylrad() { return m_cylrad; }

private:
  /** number of spheres. */
  int m_nsphere;
  /** dimension of the maze. */
  double m_dim;
  /** sphere radius. */
  double m_sphrad;
  /** cylinder (connecting the spheres) radius*/
  double m_cylrad;
};
}

#endif
