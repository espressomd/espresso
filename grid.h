#ifndef GRID_H
#define GRID_H
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

#endif
