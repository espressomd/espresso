#ifndef NSQUARE_H
#define NSQUARE_H
#include "cells.h"

Cell *nsq_position_to_cell(double pos[3]);

void nsq_topology_release();

void nsq_topology_init(CellPList *local);

void nsq_balance_particles();

#endif
