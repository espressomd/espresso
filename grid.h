#ifndef GRID_H
#define GRID_H
/** \file grid.h
 *
 *  For more information on the domain decomposition, see \ref grid.c "grid.c". 
*/

#include <tcl.h>

/**************************************************
 * exported variables
 **************************************************/

/** The number of nodes in each spatial dimension. */
extern int node_grid[3];
/** position of node in node grid */
extern int node_pos[3];
/** the six nearest neighbors of a node in the node grid. */
extern int node_neighbors[6];
/** where to fold particles that leave local box in direction i. */
extern int boundary[6];
/** whether a node is infinitely extended in that direction
    (necessary for non periodic coordinates). */
extern int extended[6];

/** Simulation box dimensions. */ 
extern double box_l[3];
/** Flags for all three dimensions wether pbc are applied (default). */ 
extern int periodic[3];
/** Dimensions of the box a single node is responsible for. */ 
extern double local_box_l[3];
/** Left (bottom, front) corner of this nodes local box. */ 
extern double my_left[3];
/** Right (top, back) corner of this nodes local box. */ 
extern double my_right[3];


/**************************************************
 * functions
 **************************************************/

/** Make sure that the node grid is set, eventually
    determine one automatically. */
void setup_node_grid();

/** return wether node grid was set. */
int node_grid_is_set();

/** node mapping: array -> node. 
 *
 * @param node   number of the node you want to know the position for.
 * @param pos[3] position of the node in node grid.        
*/
void map_node_array(int node, int pos[3]);

/** node mapping: node -> array. 
 *
 * @return       number of the node at position pos.
 * @param pos[3]  position of the node in node grid.        
*/
int map_array_node(int pos[3]);

/** map particle position to node. 
 *
 * @return       number of the node where the particle is (should be) located.
 * @param pos[3] particle position
*/
int find_node(double pos[3]);

/** datafield callback for \ref node_grid. */
int node_grid_callback(Tcl_Interp *interp, void *data);

/** datafield callback for \ref box_l. Sets the box dimensions. */
int boxl_callback(Tcl_Interp *interp, void *_data);

/** datafield callback for \ref periodic. Determines wether a coordinate is pbc (default). */
int per_callback(Tcl_Interp *interp, void *_data);

/** fill neighbor lists of node. 
 *
 * Calculates the numbers of the 6 nearest neighbors for a node and
 * stors them in \ref node_neighbors.
 *
 * @param node number of the node.  */
void calc_node_neighbors(int node);

/** Call this if the topology (grid, box dim, ...) has changed. Only for master node,
    will be communicated. */
void changed_topology();

#endif
