#ifndef GHOSTS_H 
#define GHOSTS_H 
/** \file ghosts.h
    For more information on ghosts,
    see \ref ghosts.c "ghosts.c"
 */

/************************************************
 * data types
 ************************************************/

/** Structure to hold ghost particle information. 
 *  If this changes also the functions pack_gost/unpack_ghost 
 *  have to be reweritten. */
typedef struct {
  int identity;
  int type;
  double p[3];
  double q;
} Ghost;

/************************************************
 * functions
 ************************************************/

/** initialize ghost particle structures. 
 *  
 *  <ul>
 *      <li> Init node neighbours ( \ref calc_node_neighbors ).
 *      <li> Init particle exchange buffers.
 *      <li> Init ghost exchange cell structures:
 *  </ul>      
 */
void ghost_init();

/** exchange particles between neighbouring nodes. 
    For a new verlet list setup first all local particles which have left the local box have to be send to the right node. Procedur: direction_loop: 1) 
*/
void exchange_part();
/** exchange ghost particles. */
void exchange_ghost();
/** exchange ghost particle positions. */
void update_ghost_pos();
/** exchange ghost particle forces. */
void collect_ghost_forces();
/** exit ghost structures. */
void ghost_exit();

#endif
