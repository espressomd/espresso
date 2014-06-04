/*
  Copyright (C) 2010,2012,2013 The ESPResSo project
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
/** \file layered.hpp
    The layered cellsystem. This cellsystem is a combination of a single processor n-squared method along x and y,
    and a multiprocessor domain decomposition along z. Therefore only \f$1\times 1\times N\f$ processors grids are
    allowed for this cellsystem. The implementation is pretty similar to
    \ref domain_decomposition.hpp "domain_decomposition.h".
*/
#ifndef LAYERED_H
#define LAYERED_H
#include "cells.hpp"

/** number of layers, i. e. cells, per node */
extern int n_layers, determine_n_layers;

/** height of the layers, i. e. box_l[2]/n_nodes */
extern double layer_h, layer_h_i;

/** map a position to a cell, if on this node, else returns NULL. */
Cell *layered_position_to_cell(double pos[3]);

/// free all data structure that belong to this cell system
void layered_topology_release();

/// initialize the layered cell system and sort in the particles
void layered_topology_init(CellPList *local);

/// distribute all particles such that they are in their dedicated cell
void layered_exchange_and_sort_particles(int global_flag);

/// calculate short ranged forces
void layered_calculate_ia();

/// calculate short ranged energy
void layered_calculate_energies();

/// calculate short ranged virials
void layered_calculate_virials(int v_comp);

/// calculate the minimum image vector
void layered_get_mi_vector(double res[3], double a[3], double b[3]);

#endif
