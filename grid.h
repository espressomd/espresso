#ifndef GRID_H
#define GRID_H
/** \file grid.h   Domain decomposition for parallel computing.
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:limbach@mpip-mainz.mpg.de">Hanjo</a>
 *
 *  The primary simulation box is divided into orthogonal rectangular
 *  subboxes which are assigned to the different nodes (or processes
 *  or processors or threads if you want). This grid is described in
 *  \ref node_grid. Each node has a number \ref this_node and a
 *  position \ref node_pos in that grid. Each node has also 6 nearest
 *  neighbors \ref node_neighbors which are necessary for the
 *  communication between the nodes (see also \ref ghosts.c and \ref
 *  p3m.c for more details about the communication.
 *
 *  For the 6 directions \anchor directions we have the following convention:
 *
 *  \image html directions.gif "Convention for the order of the directions"
 *  \image latex directions.eps "Convention for the order of the directions" width=6cm
 *
 *  The Figure illustrates the direction convetion used for arrays
 *  with 6 (e.g. \ref node_neighbors, \ref boundary) and 3 entries
 *  (e.g \ref node_grid, \ref box_l , \ref my_left,...).
 *  
 *  Attention: If you change anything of the simulation box dimensions
 *  you have to call \ref changed_topology.
 *
 *  For more information on the domain decomposition, see \ref grid.c "grid.c". 
*/
#include <tcl.h>

/** \name Exported Variables */
/************************************************************/
/*@{*/

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
/** Flags for all three dimensions wether pbc are applied (default). */ 
extern int periodic[3];

/** Simulation box dimensions. */ 
extern double box_l[3];
/** Dimensions of the box a single node is responsible for. */ 
extern double local_box_l[3];
/** Left (bottom, front) corner of this nodes local box. */ 
extern double my_left[3];
/** Right (top, back) corner of this nodes local box. */ 
extern double my_right[3];

/*@}*/

/** \name Exported Functions */
/************************************************************/
/*@{*/

/** Make sure that the node grid is set, eventually
    determine one automatically. */
void setup_node_grid();

/** return wether node grid was set. */
int node_grid_is_set();

/** node mapping: array -> node. 
 *
 * \param node   number of the node you want to know the position for.
 * \param pos[3] position of the node in node grid.        
*/
void map_node_array(int node, int pos[3]);

/** node mapping: node -> array. 
 *
 * \return       number of the node at position pos.
 * \param pos[3]  position of the node in node grid.        
*/
int map_array_node(int pos[3]);

/** map particle position to node. 
 *
 * \return       number of the node where the particle is (should be) located.
 * \param pos[3] particle position
*/
int find_node(double pos[3]);

/** fill neighbor lists of node. 
 *
 * Calculates the numbers of the 6 nearest neighbors for a node and
 * stors them in \ref node_neighbors.
 *
 * \param node number of the node.  */
void calc_node_neighbors(int node);

/** Call this if the topology (grid, box dim, ...) has changed. Only for master node,
    will be communicated. */
void changed_topology();

/** calculate most square 2d grid. */
void calc_2d_grid(int n, int grid[3]);

/** Calculate most cubic 3d grid. */
void calc_3d_grid(int n, int grid[3]);

/** calculate 'best' mapping between a 2d and 3d grid.
 *  This we need for the communication from 3d domain decomposition 
 *  to 2d row decomposition. 
 *  The dimensions of the 2d grid are resorted, if necessary, in a way
 *  that they are multiples of the 3d grid dimensions.
 *  \param g3d      3d grid.
 *  \param g2d      2d grid.
 *  \param mult     factors between 3d and 2d grid dimensions
 *  \return         index of the row direction [0,1,2].
*/ 
int map_3don2d_grid(int g3d[3],int g2d[3], int mult[3]);

/** datafield callback for \ref node_grid. */
int node_grid_callback(Tcl_Interp *interp, void *data);

/** datafield callback for \ref box_l. Sets the box dimensions. */
int boxl_callback(Tcl_Interp *interp, void *_data);

/** datafield callback for \ref periodic. Determines wether a coordinate is pbc (default). */
int per_callback(Tcl_Interp *interp, void *_data);

/*@}*/
#endif
