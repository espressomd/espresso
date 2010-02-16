// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2009; all rights reserved unless otherwise stated.
/** \file layered.h
    The layered cellsystem. This cellsystem is a combination of a single processor n-squared method along x and y,
    and a multiprocessor domain decomposition along z. Therefore only \f$1\times 1\times N\f$ processors grids are
    allowed for this cellsystem. The implementation is pretty similar to
    \ref domain_decomposition.h "domain_decomposition.h".
*/
#ifndef LAYERED_H
#define LAYERED_H
#include "cells.h"

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
void layered_calculate_virials();

/// calculate the minimum image vector
void layered_get_mi_vector(double res[3], double a[3], double b[3]);

#endif
