/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
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
/** \file lattice_inline.hpp */

#ifndef LATTICE_INLINE_HPP
#define LATTICE_INLINE_HPP

#include "lattice.hpp"

inline int Lattice::map_lattice_to_node(int *ind, int *grid) {

  /* determine coordinates in node_grid */
  grid[0] = (int)floor(ind[0]*this->agrid[0]*box_l_i[0]*node_grid[0]);
  grid[1] = (int)floor(ind[1]*this->agrid[1]*box_l_i[1]*node_grid[1]);
  grid[2] = (int)floor(ind[2]*this->agrid[2]*box_l_i[2]*node_grid[2]);

  //fprintf(stderr,"%d: (%d,%d,%d)\n",this_node,grid[0],grid[1],grid[2]);

  /* change from global to local lattice coordinates */
  ind[0] = ind[0] - grid[0]*this->grid[0] + this->halo_size;
  ind[1] = ind[1] - grid[1]*this->grid[1] + this->halo_size;
  ind[2] = ind[2] - grid[2]*this->grid[2] + this->halo_size;

  /* return linear index into node array */
  return map_array_node(grid);
}

#endif // LATTICE_INLINE_HPP
