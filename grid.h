/************************************************/
/*******************  GRID.H  *******************/
/************************************************/
#ifndef GRID_H
#define GRID_H

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "global.h"

/** try to determine the processor grid and communicate it */
int setup_processor_grid();

/** return wether grid was set */
int processor_grid_is_set();

/** node mapping array -> node */
void map_node_array(int node, int *a, int *b, int *c);

/** node mapping node -> array */
int map_array_node(int a, int b, int c);

/** map position to node */
int find_node(double pos[3]);

/** datafield callback for procgrid */
int pgrid_callback(Tcl_Interp *interp, void *data);

#endif
