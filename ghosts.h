// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2003; all rights reserved unless otherwise stated.
#ifndef GHOSTS_H 
#define GHOSTS_H 
/** \file ghosts.h    Ghost particles and particle exchange.
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:limbach@mpip-mainz.mpg.de">Hanjo</a>
 *
 *  In this file you find everything concerning the exchange of
 *  particle data (particles, ghosts, positions and forces) for short
 *  range interactions between the spacial domains of neighbouring
 *  nodes.
 *
 *  All structures are initialized by \ref ghost_pre_init and \ref
 *  ghost_init.
 *
 *  Exchanging particles that move from one cell to another including
 *  crossing of node domain boundaries is done by \ref
 *  exchange_and_sort_part. 
 *
 *  Setup of the ghost particle frame is done by \ref exchange_ghost.
 *
 *  During integration the ghost positions and forces are exchanged by
 *  \ref update_ghost_pos and \ref collect_ghost_forces.
 *
 *  Communication \anchor communication is done only between
 *  neighboring nodes. The communication scheme can be seen in the
 *  figure below.
 *
 *  \image html ghost_communication.gif "Scheme of ghost/particle communication"
 *  \image latex ghost_communication.eps "Scheme of ghost/particle communication" width=8cm 
 *
 *  To reduce the number of communications per node from 26 (number of
 *  neighbor nodes in 3 dimensions) to 6 the order of the
 *  communication is important:
 *  <ol> 
 *      <li> x-direction: left and right
 *      <li> y-direction: forth and back 
 *      <li> z-direction: down and up 
 *  </ol>
 *  In this way also edges and corners are communicated
 *  (See also \ref directions for our conventions).
 *
 *  \warning \b Since the ghost particle structures make use of the
 *  linked cell structure, \ref ghost_init has to be called after \ref
 *  cells_re_init 
 *
 *  For more information on ghosts,
 *  see \ref ghosts.c "ghosts.c"
 */
#include <mpi.h>
#include "particle_data.h"

/** \name Defines */
/************************************************************/
/*@{*/

#define GHOST_SEND 0
#define GHOST_RECV 1
#define GHOST_BCST 2
#define GHOST_RDCE 3
#define GHOST_JOBMASK  15
#define GHOST_PREFETCH 16

/*@}*/

/** \name Transfer data classes */
/************************************************************/
/*@{*/

#define GHOSTTRANS_PROPRTS  1
#define GHOSTTRANS_POSITION 2
#define GHOSTTRANS_MOMENTUM 4
#define GHOSTTRANS_FORCE    8
#define GHOSTTRANS_PARTNUM  16

/*@}*/

/** \name Data Types */
/************************************************************/
/*@{*/

typedef struct {

  /** Communication type. */
  int type;
  /** Node to communicate with (to use with GHOST_SEND, GHOST_RECV). */
  int node;
  /** MPI communicator handle (to use with GHOST_BCST, GHOST_GATH, GHOST_RDCE). */
  MPI_Comm mpi_comm;

  /** Number of particle lists to communicate. */
  int n_part_lists;
  /** Pointer array to particle lists to communicate. */
  ParticleList **part_lists;

} GhostCommunication;

/** Properties for a ghost communication. A ghost communication is defined */
typedef struct {

  /** Particle data parts to transfer */
  int data_parts;

  /** number of communication steps. */
  int num;

  /** List of ghost communications. */
  GhostCommunication *comm;

} GhostCommunicator;

/*@}*/

/** \name Exported Functions */
/************************************************************/
/*@{*/

/** Initialize a communicator. */
void prepare_comm(GhostCommunicator *comm, int data_parts, int num);

/** Free a communicator. */
void free_comm(GhostCommunicator *comm);

/** Initialize some arrays connected to the ghosts modul. 
 */
void ghost_pre_init();

/** Initialize ghost particle structures. 
 *  
 *  <ol>
 *      <li> Calculate number of cells to send/recv for each direction.
 *      <li> Allocate \ref #send_cells, \ref #recv_cells, 
 *           \ref #n_send_ghosts, \ref #n_recv_ghosts.
 *      <li> Fill \ref #send_cells and \ref #recv_cells with cell indices.
 *  </ol>
 *
 *   
 *  <ul>
 *      <li> \ref send_cells contains a list of tcells to send. 
 *      <li> \ref recv_cells contains a list of cells to receive. 
 *      <li> \ref n_send_ghosts contains the number of particles 
 *           in each send cell. 
 *      <li> \ref n_recv_ghosts contains the number of particles 
 *           in each recv cell. 
 *  </ul>
 *
 * (see \ref exchange_ghost for more details). 
*/
void ghost_init();

/** exchange particles between neighbouring nodes and sort particles in \ref #cells. 
 *
 *  Tasks:
 *  <ul> 
 *      <li> Particles which have left their cell have to be moved to their new cell.
 *      <li> Particles which have left the local box have to be send to the
 *           corresponding new node befor e.g. a new verlet list can be built.
 *      <li> Fold particle positions if the simulation box boundary has been crossed.
 *  </ul>
 * 
 *  Procedur:
 *  <ol>
 *      <li> direction loop d = 0,1,2 (x, y, z), lr = 0,1 (left, right) 
 *           (see \ref communication for details):
 *           <ul>
 *               <li> Move particles to send buffers;
 *                    via \ref move_to_p_buf to \ref p_send_buf and \ref b_send_buf.
 *               <li> Send/Receive Particles: \ref send_particles.
 *               <li> Append reveived particles to particle array of the corresponding cell:
 *                    \ref p_recv_buf and b_recv_buf via \ref append_particles to 
 *                    \ref Cell::pList .
 *           </ul>
 *           if \ref node_grid[d] == 1 (single node case for that
 *           direction) nothing is done here. \\
 *           if d=2 and lr=1 (last direction) then move particles to new cell 
 *           if necessary via \ref move_unindexed_particle.
 *      <li> Update local particle array via \ref update_local_particles.  
 * </ol> */
void exchange_and_sort_part();

/** exchange ghost particles. 
 *
 *  Procedure:
 *  <ul>
 *      <li> direction loop d = 0,1,2,3,4,5 (see \ref communication for details):
 *           <ul> 
 *               <li> Loop through \ref send_cells for this direction. 
 *               <li> Store number of particles in cells in \ref n_send_ghosts.
 *               <li> Realloc \ref g_send_buf.
 *               <li> Pack particles of that cells into \ref g_send_buf 
 *                    and fold ghost coordinates.
 *               <li> Fold ghost coordinates if necessary (\ref #boundary).
 *               <li> send/receive ghost particles (\ref send_ghosts)
 *               <li> Loop through \ref recv_cells for this direction. 
 *               <li> Store received ghosts: \ref g_recv_buf 
 *                    to \ref Cell::pList.
 *               <li> update \ref local_particles if necesarry (real particle priority!).
 *           </ul>
 *      <li> realloc buffers \ref send_buf and \ref recv_buf for 
 *           \ref update_ghost_pos and \ref collect_ghost_forces
 *  </ul>
 *
 * All buffers are resized dynamically.
 *
 * In the example (picture below) it is shown which cells Node 2
 * receives from Node 1 from direction 0 (left) and which cells Node 2
 * sends to Node 3 in direction 2 (down).
 *
 * \image html ghost_cells.gif "ghost exchange: cell structures"
 * \image latex ghost_cells.eps "ghost exchange: cell structures" \width=14cm
 *
*/
void exchange_ghost();

/** remove the ghost entries from \ref local_particles. */
void invalidate_ghosts();

/** exchange ghost particle positions. 
 *  
 *  Using the structures build up in \ref ghost_init and \ref
 *  exchange_ghost this routine updates the positions of the ghost
 *  particles.
 *
 *  It loops through the \ref directions, \ref send_cells and the
 *  particles therein and communicates their positions between
 *  neighbouring processors. Cells to send and receive (\ref
 *  send_cells, \ref recv_cells), as well as the number of particles
 *  in each cell to send/receive (\ref n_send_ghosts, \ref
 *  n_recv_ghosts) are allready known from \ref exchange_ghost.
 *
*/
void update_ghost_pos();

/** Collect ghost particle forces. 
 *
 *  Same procedure as in \ref update_ghost_pos except that the order
 *  of the directions is reversed and that forces are communicated and
 *  added up.  
*/
void collect_ghost_forces();

/*@}*/

#endif
