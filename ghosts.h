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

 *      <li> Init ghost exchange cell structures: The \ref ghost_exchange
 *           makes use of the linked cell structure (\ref cells) and hence
 *           needs more complicated structures for bookkeeping:
 *           <ul>
 *               <li> \ref send_cells is a list of all cells to send. 
 *               <li> \ref recv_cells is a list of all cells to receive. 
 *               <li> \ref n_send_ghosts contains the number of particles 
 *                    in each send cell. 
 *               <li> \ref n_recv_ghosts contains the number of particles 
 *                    in each recv cell. 
 *           </ul>
 *           All lists have the length \ref ntot_send_cells and are
 *           divided into blocks of cells send to/received from
 *           different directions. The block length is stored in \ref
 *           n_send_cells[d] , \ref n_recv_cells[d] and the beginning
 *           of each block is stored in \ref cell_start[d] . This is
 *           exemplary shown in the figure \anchor ghost_cells below.
 *
 * \image html ghost_cells.gif "ghost exchange: cell structures"
 * \image latex ghost_cells.eps "ghost exchange: cell structures" \width=15cm
 *
 *           In the example it is shown which cells Node 2 receives
 *           from Node 1 from direction 0 (left) and which cells Node
 *           2 sends to Node 3 in direction 2 (down).
 *           
 *           (see \ref exchange_ghost for more details).
 * </ul> 
*/
void ghost_init();

/** exchange particles between neighbouring nodes. 
 *
 *  Particles which have left the local box have to be send to the
 *  right node befor e.g. a new verlet list can be built.
 * 
 *  Procedur:
 *  <ol>
 *      <li> \ref fold_particle
 *      <li> direction loop d = 0,1,2 (see \ref communication for details):
 *          <ul>
 *              <li> Move particles to send buffers: \ref particles
 *                   via \ref move_to_p_buf to \ref p_send_buf and \ref b_send_buf.
 *              <li> Send/Receive Particles: \ref send_particles.
 *              <li> Move reveived particles to local particle field:
 *                   \ref p_recv_buf and b_recv_buf via \ref append_particles to 
 *                   \ref particles .
 *          </ul>
 *          if \ref node_grid[d] == 1 (single node case for that
 *          direction) nothing is done here.  
 * </ol> */
void exchange_part();

/** exchange ghost particles. 
 *
 *  Procedure:
 *  <ul>

 *      <li> Remove all ghosts from local particle array (\ref
 *           particles): \ref n_ghosts=0.

 *      <li> direction loop d = 0,1,2,3,4,5 (see \ref communication for details):
 *           <ul> 
 *               <li> Loop through \ref send_cells for this direction. 
 *               <li> Store number of particles in that cell: 
 *                    \ref n_send_ghosts[ind] = 
 *                    \ref cells[\ref send_cells[ind]].n_particles .
 *               <li> Pack particles of that cells into \ref g_send_buf via 
 *                    \ref pack_ghost.
 *               <li> Store number of ghosts to send in that direction in 
 *                    \ref ghost_send_size[d].
 *               <li> Fold ghost coordinates if necessary (\ref boundary).
 *               <li> send/receive ghost particles (\ref send_ghosts)
 *               <li> Loop through \ref recv_cells for this direction. 
 *               <li> Store received ghosts: \ref g_recv_buf via \ref unpack_ghost
 *                    to \ref particles (see also \ref particle_array).
 *           </ul>
 *  </ul>
 *  All buffers are resized dynamically.
*/
void exchange_ghost();

/** exchange ghost particle positions. 
 *  
 *  Using the structures build up in \ref ghost_init and \ref
 *  exchange_ghost thsi routine updates the positions of the ghost
 *  particles.
 *
 *  It loops through the \ref directions, \ref send_cells and the
 *  particles therein and communicates them between neighbouring
 *  processors. Cells to send and receive (\ref send_cells, \ref
 *  recv_cells), as well as the number of particles in each cell to
 *  send/receive (\ref send_ghosts, \ref recv_ghosts) are allready
 *  known from \ref exchange_ghosts.
 *
*/
void update_ghost_pos();
/** exchange ghost particle forces. 
 *
 *  Same procedure as in \ref update_ghost_pos except that the order
 *  of the directions is reversed and that forces are communicated and
 *  added up.  
*/
void collect_ghost_forces();
/** exit ghost structures. */
void ghost_exit();

#endif
