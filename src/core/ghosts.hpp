/*
  Copyright (C) 2010,2011,2012,2013 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
    Max-Planck-Institute for Polymer Research, Theory Group
  
  This file is part of ESPResSo.
  
  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/
#ifndef _GHOSTS_H 
#define _GHOSTS_H 
/** \file ghosts.hpp    Ghost particles and particle exchange.
 
In this file you find everything concerning the exchange of
particle data (particles, ghosts, positions and forces) for short
range interactions between the spacial domains of neighbouring
nodes.

<h2> How does this work </h2>
The ghost communication transfers data from cells on one node to cells on another node during
the integration process. Note that data can only be transferred from one cell to one other, and
the contents of the other cell will be overwritten. This communication is invoked by the integrator,
at four different times in an integration step. These steps are reflected in \ref cell_structure structure,
since they are treated differently by different cell systems:
<ul>
<li> ghost_cells is used to transfer the cell sizes, i. e. makes sure that for all later transfers
      the target cell has the same size as the source cell.
<li> exchange_ghosts sets up all information on the ghosts that is necessary. Normally transfers the
      (shifted) position and particle properties.
<li> update_ghosts is used to update the particle properties if no particles have been moved between cells.
      Therefore only the positions are transferred, otherwise this normally looks pretty much like the
      exchange_ghosts.
<li> collect_ghost_forces finally is used to transfer back forces that were exerted on a ghost particle. They
      are simply added again to the force of the real particle. The communication process is therefore inverted
      with respect to exchange_ghosts and update_ghosts, i. e. sending is replaced by receiving and the other way
      round.
</ul>
The particle data that has to be transferred, and especially from where to where, heavily depends on the cell system.
In Espresso, this is abstracted in form of ghost communicators (\ref GhostCommunicator) and ghost communications
(\ref GhostCommunication). The ghost communicators represent the four communications above and consist of the data to
transfer (which is determined by their type) and a list of ghost communications.
The data types are described by the particle data classes
<ul>
<li> GHOSTTRANS_PROPRTS transfers the \ref ParticleProperties
<li> GHOSTTRANS_POSITION transfers the \ref ParticlePosition
<li> GHOSTTRANS_POSSHFTD is no transfer class, but marks that an additional vector, the shift, has to be
      added to the positions. Can only be used together with GHOSTTRANS_POSITION.
<li> GHOSTTRANS_MOMENTUM transfers the \ref ParticleMomentum
<li> GHOSTTRANS_FORCE transfers the \ref ParticleForce
<li> GHOSTTRANS_PARTNUM transfers the cell sizes
</ul>
Each ghost communication describes a single communication of the local with another node (or all other nodes). The data
transferred can be any number of cells, there are five communication types:
<ul>
<li> GHOST_SEND sends data to one other node
<li> GHOST_RECV recvs data from one other node. In the case of forces, they are added up, everything else is overwritten
<li> GHOST_BCST sends data to all nodes
<li> GHOST_RDCE recvs data from all nodes.  In the case of forces, they are added up, everything else is overwritten
<li> GHOST_LOCL transfer data from a local cell to another local cell. In this case, the first half of the cells are the sending cells,
      the other half are the corresponding receivers.
</ul>
Note that for the first four communications you have to make sure that when one node sends to another, that the sender has GHOST_SEND and the
receiver GHOST_RECV. In the case of GHOST_BCST rsp. GHOST_RDCE, all nodes have to have the same communication type and the same master
sender/receiver (just like the MPI commands).

A special topic are GHOST_PREFETCH and GHOST_PSTSTORE. For example, if all nodes broadcast to the other, the naive implementation will be that
n_nodes times a GHOST_BCST is done with different master nodes. But this means that each time n_nodes-1 nodes wait for the master to construct
its send buffer. Therefore there is the prefetch flag which can be set on a pair of recv/send operations. If the ghost communication reaches a
recv operation with prefetch, the next send operation (which must have the prefetch set!!) is searched and the send buffer already created.
When sending, this precreated send buffer is used. In the scenario above, all nodes create the send buffers simultaneously in the first
communication step, thereby reducing the latency a little bit. The pststore is similar and postpones the write back of received data until a
send operation (with a precreated send buffer) is finished.

The ghost communicators are created in the init routines of the cell systems, therefore have a look at \ref dd_topology_init or
\ref nsq_topology_init for further details.
*/
#include <mpi.h>
#include "particle_data.hpp"

/** \name Transfer types, for \ref GhostCommunicator::type */
/************************************************************/
/*@{*/

/// send to a single node
#define GHOST_SEND 0
/// recv from a single node
#define GHOST_RECV 1
/// broadcast, the node entry gives the sender
#define GHOST_BCST 2
/// reduce, the node entry gives the receiver
#define GHOST_RDCE 3
/// transfer data from cell to cell on this node
#define GHOST_LOCL 4

/// mask to the job area of the transfer type
#define GHOST_JOBMASK  15
/// additional flag for prefetching
#define GHOST_PREFETCH 16
/// additional flag for poststoring
#define GHOST_PSTSTORE 32
/*@}*/


/** \name Transfer data classes, for \ref GhostCommunication::type */
/************************************************************/
/*@{*/
/// transfer \ref ParticleProperties
#define GHOSTTRANS_PROPRTS  1
/// transfer \ref ParticlePosition
#define GHOSTTRANS_POSITION 2
/** flag for \ref GHOSTTRANS_POSITION, shift the positions by \ref GhostCommunication::shift.
    Must be or'd together with \ref GHOSTTRANS_POSITION */
#define GHOSTTRANS_POSSHFTD 4
/// transfer \ref ParticleMomentum
#define GHOSTTRANS_MOMENTUM 8
/// transfer \ref ParticleForce
#define GHOSTTRANS_FORCE    16

#ifdef LB
/// transfer \ref ParticleLatticeCoupling
#define GHOSTTRANS_COUPLING 32
#endif

/// resize the receiver particle arrays to the size of the senders
#define GHOSTTRANS_PARTNUM  64
/*@}*/

/** \name Data Types */
/************************************************************/
/*@{*/

typedef struct {

  /** Communication type. */
  int type;
  /** Node to communicate with (to use with all MPI operations). */
  int node;
  /** MPI communicator handle (to use with GHOST_BCST, GHOST_GATH, GHOST_RDCE). */
  MPI_Comm mpi_comm;

  /** Number of particle lists to communicate. */
  int n_part_lists;
  /** Pointer array to particle lists to communicate. */
  ParticleList **part_lists;

  /** if \ref GhostCommunicator::data_parts has \ref GHOSTTRANS_POSSHFTD, then this is the shift vector.
      Normally this a integer multiple of the box length. The shift is done on the sender side */
  double shift[3];
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

/** Initialize ghosts. */
void ghost_init();

/** do a ghost communication */
void ghost_communicator(GhostCommunicator *gc);

/** Go through \ref ghost_cells and remove the ghost entries from \ref
    local_particles. Part of \ref dd_exchange_and_sort_particles.*/
void invalidate_ghosts();

/* TODO: This function is not used anywhere. To be removed?  */
#ifdef GHOST_FLAG
inline int ifParticleIsGhost(Particle *p){
   return p->l.ghost;
}
#endif

/*@}*/

#endif
