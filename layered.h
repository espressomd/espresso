#ifndef LAYERED_H
#define LAYERED_H
#include "cells.h"

#define LAYERED_FULL_EXCHANGE 0
#define LAYERED_NBOR_EXCHANGE 1
extern int n_layers;

Cell *layered_position_to_cell(double pos[3]);

void layered_topology_release();

void layered_topology_init(CellPList *local);

void layered_exchange_and_sort_particles(int global_flag);

void layered_calculate_ia();

void layered_calculate_energies();

void layered_calculate_virials();

#endif
