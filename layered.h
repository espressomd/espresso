#ifndef LAYERED_H
#define LAYERED_H
#include "cells.h"

extern int n_layers, determine_n_layers;

extern double layer_h, layer_h_i;

Cell *layered_position_to_cell(double pos[3]);

void layered_topology_release();

void layered_topology_init(CellPList *local);

void layered_exchange_and_sort_particles(int global_flag);

void layered_calculate_ia();

void layered_calculate_energies();

void layered_calculate_virials();

#endif
