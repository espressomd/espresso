/************************************************/
/*******************  GRID.H  *******************/
/************************************************/
#ifndef GRID_H
#define GRID_H
#include <tcl.h>
/**************************************************
 * exported variables
 **************************************************/

/** The number of nodes in each spatial dimension. */
extern int processor_grid[3];
extern int pe_pos[3];
extern int neighbors[6];
extern double boundary[6];

/** Simulation box dimensions. */ 
extern double box_l[3];
/** Dimensions of the box a single node is responsible for. */ 
extern double local_box_l[3];
/** Left top corner of this nodes local box. */ 
extern double my_left[3];
/** Right bottom corner of this nodes local box. */ 
extern double my_right[3];


/**************************************************
 * functions
 **************************************************/

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

/** fill neighbor lists of node. */
void calc_neighbors(int node);

/** call if topology (grid, box dim, ...) changed */
void changed_topology();

#endif
